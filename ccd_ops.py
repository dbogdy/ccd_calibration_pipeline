#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ccdproc-based pixel operations shared by the pipeline scripts (v3).

All image arithmetic goes through ccdproc / astropy.nddata, following the
Astropy "CCD Data Reduction Guide":

  - reading      : astropy.nddata.CCDData (unit from BUNIT, 'adu' fallback)
  - overscan     : ccdproc.subtract_overscan (BIASSEC) + trim_image (TRIMSEC)
  - combination  : ccdproc.combine, method='average', sigma clipping with
                   np.ma.median / astropy.stats.mad_std and a memory limit
                   (the guide's recommended master-frame recipe)
  - calibration  : ccdproc.subtract_bias / subtract_dark / flat_correct
                   (used directly in the step scripts)

The per-frame statistics / quality-control helpers work on plain numpy
arrays (they are diagnostics, not calibration math).
"""

import tempfile
import warnings
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import astropy.units as u
import ccdproc as ccdp
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import mad_std, sigma_clipped_stats
from astropy.wcs import FITSFixedWarning

from calib_common import MEM_BUDGET_MB, env_choice, env_float, env_int

# CCDData.read parses the WCS, and astropy's datfix "repairs" headers that
# have DATE-OBS but no MJD-OBS - a harmless completion that would otherwise
# print one FITSFixedWarning per file read. Irrelevant for EL imaging.
warnings.filterwarnings("ignore", category=FITSFixedWarning)

DEFAULT_UNIT = "adu"

# All thresholds can be overridden from config.ini via run_pipeline.py
# (the config sections map to CALIB_* environment variables); standalone
# runs use these defaults.

# Sigma-clipping thresholds for master combination, per the Astropy CCD
# guide (5 sigma around the median, MAD as the robust deviation).
SIGMA_CLIP_LOW = env_float("CALIB_SIGMA_CLIP_LOW", 5.0)
SIGMA_CLIP_HIGH = env_float("CALIB_SIGMA_CLIP_HIGH", 5.0)

# How the master builders combine their frames. 'average' (with the sigma
# clipping above) is the guide's recommendation - lower noise than the
# median at equal robustness; 'median' is available for special cases.
MASTER_COMBINE_METHOD = env_choice("CALIB_COMBINE_METHOD", "average",
                                   ("average", "median"))

# Groups smaller than this do not produce a master at all. Distinct from
# MIN_FRAMES_FOR_QC, which only decides whether outliers may be rejected:
# a "master" combined from one or two frames is essentially a raw frame,
# keeping its full read noise and gaining nothing from sigma clipping.
MIN_FRAMES_FOR_MASTER = env_int("CALIB_MIN_FRAMES_FOR_MASTER", 5)

# QC thresholds (level-based QC, used for bias and dark frames)
MIN_FRAMES_FOR_QC = env_int("CALIB_MIN_FRAMES_FOR_QC", 5)
# below this, only warn - never reject
MEDIAN_SIGMA = env_float("CALIB_MEDIAN_SIGMA", 5.0)
# reject if |median - group| > MEDIAN_SIGMA * MAD + floor
MEDIAN_FLOOR_ADU = env_float("CALIB_MEDIAN_FLOOR_ADU", 1.0)
# absolute floor so tiny MADs don't cause false rejects
STD_RATIO_MAX = env_float("CALIB_STD_RATIO_MAX", 2.0)
# reject if frame noise > 2x group, or < 1/2x

MAX_RN_PAIRS = env_int("CALIB_MAX_RN_PAIRS", 8)
# frame pairs used for the pair-difference noise

# Per-frame statistics feed the quality control and the bias-level screen
# only - never the pixel data that ends up in a master - so they may be
# computed on a subsampled grid ([::N, ::N]). N=1 restores full-array stats.
STATS_SUBSAMPLE = max(1, env_int("CALIB_STATS_SUBSAMPLE", 4))

# Frame reads dominate the runtime (especially through WSL's /mnt drvfs
# layer), so scanning and per-frame reduction run on a small thread pool.
# 1 = fully sequential.
IO_WORKERS = max(1, env_int("CALIB_IO_WORKERS", 4))

# Optional light-frame cleaning (used by light_calibration.py). Disabled by
# default; turned on per run via config.ini [calibration]
# fix_hot_pixels / fix_cosmic_rays.
HOT_PIXEL_SIGMA = env_float("CALIB_HOT_PIXEL_SIGMA", 5.0)
# hot-pixel threshold: clipped_median + sigma * clipped_std of the master dark
CR_SIGCLIP = env_float("CALIB_CR_SIGCLIP", 4.5)
# L.A.Cosmic detection threshold [sigma above background]
CR_OBJLIM = env_float("CALIB_CR_OBJLIM", 5.0)
# L.A.Cosmic min contrast between a cosmic ray and an underlying source
CR_READNOISE = env_float("CALIB_CR_READNOISE", 8.0)
# fallback read noise [e-] used when the frame header has no RDNOISE


def read_ccd(path):
    """
    Read a FITS file as CCDData. The unit comes from BUNIT; when the header
    has no (valid) unit, raw camera counts are assumed ('adu').

    The unit is decided from the header *before* reading, rather than by
    letting CCDData.read fail on a missing BUNIT and retrying with a unit.
    That retry is not safe to run from several threads: astropy's reader can
    then return a CCDData with unit=None instead of raising ValueError, so
    the fallback never fires and the unitless frame only blows up later, in
    the first ccdproc subtraction, as a confusing UnitConversionError
    ("'adu' and '' (dimensionless) are not convertible").
    """
    bunit = None
    try:
        raw_unit = fits.getheader(path, ext=0).get("BUNIT")
        if raw_unit and str(raw_unit).strip():
            u.Unit(str(raw_unit).strip())      # validate before trusting it
            bunit = str(raw_unit).strip()
    except Exception:
        bunit = None
    # A valid BUNIT means CCDData.read finds a unit on its own and never
    # takes the failing branch; otherwise the unit is supplied explicitly.
    if bunit:
        return CCDData.read(path)
    return CCDData.read(path, unit=DEFAULT_UNIT)


def apply_overscan_trim(ccd):
    """
    Standard first step from the CCD guide: subtract the overscan level
    (BIASSEC, per-line median) and trim to the science region (TRIMSEC).
    Both are skipped silently when the keywords are absent. Returns
    (ccd, overscan_used).
    """
    used = False
    biassec = ccd.header.get("BIASSEC")
    if biassec:
        try:
            ccd = ccdp.subtract_overscan(ccd, fits_section=str(biassec),
                                         overscan_axis=None, median=True)
            used = True
        except Exception as exc:
            print(f"  WARNING: overscan (BIASSEC '{biassec}') failed: {exc}")
    trimsec = ccd.header.get("TRIMSEC")
    if trimsec:
        try:
            ccd = ccdp.trim_image(ccd, fits_section=str(trimsec))
        except Exception as exc:
            print(f"  WARNING: trim (TRIMSEC '{trimsec}') failed: {exc}")
    return ccd, used


def compact(ccd):
    """
    Prepare a CCDData for writing as a simple single-HDU float32 file:
    cast the data and drop mask/uncertainty (the rest of the pipeline reads
    plain ext=0 images). Masked pixels, if any, are set to NaN first so
    they stay visible instead of carrying arbitrary filled values.
    """
    data = np.asarray(ccd.data, dtype=np.float32)
    if ccd.mask is not None and ccd.mask.any():
        data = data.copy()
        data[ccd.mask] = np.nan
    ccd.data = data
    ccd.mask = None
    ccd.uncertainty = None
    return ccd


def warn_if_nonfinite(ccd):
    """
    Warn (on screen and in the header HISTORY) when a combined master
    contains non-finite pixels - a NaN in a master propagates into every
    frame calibrated with it. The data is NOT modified.
    """
    n_bad = int(np.count_nonzero(~np.isfinite(ccd.data)))
    if n_bad:
        warn = (f"master contains {n_bad} non-finite pixel(s) (NaN/inf) - "
                "check input frames")
        print(f"  WARNING: {warn}")
        ccd.header.add_history(f"WARNING: {warn}")
    return n_bad


def inv_median(a):
    """Per-frame scale function for flats (CCD guide): 1 / median."""
    return 1.0 / np.median(a)


def combine_frames(paths, method=None, sigma_clip=True, scale=None,
                   sigma_low=SIGMA_CLIP_LOW, sigma_high=SIGMA_CLIP_HIGH):
    """
    Combine FITS files into one master with ccdproc.combine, using the CCD
    guide's recipe: average (default; MASTER_COMBINE_METHOD, configurable)
    with sigma clipping around the median (MAD as deviation) and a memory
    limit so large stacks are processed in chunks. Pass method explicitly
    ('average'/'median'/'sum') to bypass the configured default.

    Returns a compacted CCDData (float32, no mask/uncertainty) with the
    method recorded in COMBMETH.
    """
    if method is None:
        method = MASTER_COMBINE_METHOD
    paths = [str(p) for p in paths]
    kwargs = {}
    # Raw camera files usually have no BUNIT - tell ccdproc what they are.
    if "BUNIT" not in fits.getheader(paths[0], ext=0):
        kwargs["unit"] = DEFAULT_UNIT
    master = ccdp.combine(
        paths,
        method=method,
        scale=scale,
        sigma_clip=sigma_clip,
        sigma_clip_low_thresh=sigma_low,
        sigma_clip_high_thresh=sigma_high,
        sigma_clip_func=np.ma.median,
        sigma_clip_dev_func=mad_std,
        mem_limit=MEM_BUDGET_MB * 1024 * 1024,
        **kwargs,
    )
    master.meta["COMBMETH"] = (method, "ccdproc combine method")
    return compact(master)


class ReducedTempDir:
    """
    Context manager for the guide's "write the per-frame reduced images,
    then combine the files" workflow, without keeping the intermediates:
    a temporary folder (on the same drive as the masters) that is removed
    when the group is done.
    """

    def __init__(self, base_dir: Path, label: str):
        base_dir.mkdir(parents=True, exist_ok=True)
        self._tmp = tempfile.TemporaryDirectory(
            prefix=f"reduced_{label}_", dir=str(base_dir))
        self.path = Path(self._tmp.name)

    def write(self, ccd, name: str) -> Path:
        out = self.path / name
        compact(ccd).write(out, overwrite=True)
        return out

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self._tmp.cleanup()
        return False


# ---------------------------------------------------------------------------
# Per-frame statistics and quality control (diagnostics, plain numpy)
# ---------------------------------------------------------------------------

def map_frames(func, items, workers=None):
    """
    Apply func to every item, on a thread pool when it pays off. Results keep
    the input order; frame reads release the GIL, so threads help even though
    the work is partly numeric. Used for the per-frame scans and for the
    per-frame reduction loops in the master builders.
    """
    items = list(items)
    n_workers = IO_WORKERS if workers is None else max(1, int(workers))
    if n_workers <= 1 or len(items) < 2:
        return [func(it) for it in items]
    with ThreadPoolExecutor(max_workers=min(n_workers, len(items))) as pool:
        return list(pool.map(func, items))


def _scan_one(path, step):
    """Per-frame statistics for one file; errors are returned, not raised."""
    try:
        data = fits.getdata(path, ext=0, memmap=False)
    except Exception as exc:
        return {"file": Path(path), "error": str(exc)}
    shape = data.shape
    sample = data[::step, ::step] if step > 1 else data
    sample = np.asarray(sample, dtype=np.float32)
    return {
        "file": Path(path),
        "shape": shape,
        "median": float(np.median(sample)),
        "std": float(mad_std(sample)),
    }


def scan_frames(paths, subsample=None, workers=None):
    """
    Collect per-frame statistics (median / noise). Statistics are measured on
    a [::N, ::N] grid (STATS_SUBSAMPLE, config.ini [master_settings]
    stats_subsample) because they only drive quality control, never the
    combined pixel data. Returns a list of dicts: file, shape, median, std.
    """
    step = STATS_SUBSAMPLE if subsample is None else max(1, int(subsample))
    results = map_frames(lambda p: _scan_one(p, step), paths, workers)
    frames = []
    for res in results:
        if "error" in res:
            print(f"  WARNING: cannot read data of {res['file'].name}: "
                  f"{res['error']}")
            continue
        frames.append(res)
    return frames


def _shape_one(path):
    """Frame shape straight from the header - no pixel data is read."""
    try:
        header = fits.getheader(path, ext=0)
        return {"file": Path(path),
                "shape": (int(header["NAXIS2"]), int(header["NAXIS1"]))}
    except Exception as exc:
        return {"file": Path(path), "error": str(exc)}


def scan_shapes(paths, workers=None):
    """
    Like scan_frames but reads only the headers, for callers that need the
    frame shape and nothing else (the master flat builder measures its
    statistics after calibration, so scanning the raw pixels is wasted I/O).
    Returns a list of dicts: file, shape.
    """
    results = map_frames(_shape_one, paths, workers)
    frames = []
    for res in results:
        if "error" in res:
            print(f"  WARNING: cannot read header of {res['file'].name}: "
                  f"{res['error']}")
            continue
        frames.append(res)
    return frames


def enough_frames(count, files=None, stage="", minimum=None):
    """
    True when a group still has enough frames to be worth combining.

    A group below the threshold is almost always a metadata typo - one stray
    GAIN / XBINNING / SET-TEMP value splitting a handful of frames off the
    group they belong to - so the offending files are listed. Building a
    master from one or two frames would write a file that is essentially a
    raw frame, with its full read noise and nothing for sigma clipping to
    reject, under a name the master lookup would happily select later.
    """
    limit = MIN_FRAMES_FOR_MASTER if minimum is None else minimum
    if count >= limit:
        return True
    where = f" {stage}" if stage else ""
    print(f"  SKIPPED      : only {count} frame(s){where}, below "
          f"min_frames_for_master ({limit}) - no master written")
    if files:
        files = list(files)
        for f in files[:5]:
            print(f"  {'':<13}  {Path(f).name}")
        if len(files) > 5:
            print(f"  {'':<13}  ... and {len(files) - 5} more")
        print(f"  {'':<13}  check GAIN / XBINNING / SET-TEMP on these - a "
              f"group this small is usually a header typo")
    return False


def filter_majority_shape(frames):
    """
    All frames in a group must share one shape; drop the odd ones out
    (with a warning). Returns (frames, shape).
    """
    shapes = [fr["shape"] for fr in frames]
    shape = max(set(shapes), key=shapes.count)
    for fr in [fr for fr in frames if fr["shape"] != shape]:
        print(f"  WARNING: {fr['file'].name} has shape {fr['shape']}, "
              f"group uses {shape} - excluded")
    return [fr for fr in frames if fr["shape"] == shape], shape


def print_group_stats(frames, label="Frame"):
    """
    One-line-per-quantity summary of the group (levels and noise), instead of
    a table with one row per frame.
    """
    if not frames:
        return
    meds = np.array([fr["median"] for fr in frames], dtype=float)
    stds = np.array([fr["std"] for fr in frames], dtype=float)
    print(f"  {label + ' levels':<13}: median {np.median(meds):.2f} ADU  "
          f"(range {meds.min():.2f} .. {meds.max():.2f}, "
          f"MAD {mad_std(meds):.2f})")
    print(f"  {label + ' noise':<13}: median {np.median(stds):.3f} ADU  "
          f"(range {stds.min():.3f} .. {stds.max():.3f})")


def print_qc_summary(kept, rejected, info, max_listed=3):
    """
    Compact quality-control report: the counts and the thresholds actually
    used, then at most max_listed rejected frames with the rest aggregated
    by reason.
    """
    limits = []
    if info.get("median_tol") is not None:                 # bias / dark QC
        limits.append(f"level tol +-{info['median_tol']:.2f}")
    if info.get("std_lo") is not None:
        limits.append(f"noise [{info['std_lo']:.2f}, {info['std_hi']:.2f}]")
    if info.get("low_limit") is not None:                  # flat QC
        limits.append(f"median [{info['low_limit']:.0f}, "
                      f"{info['sat_limit']:.0f}] ADU")
    detail = f"   ({', '.join(limits)})" if limits else ""
    print(f"  {'QC':<13}: {kept} kept, {len(rejected)} rejected{detail}")

    if info.get("safety_net"):
        print(f"  {'':<13}  WARNING: QC would have rejected most of the group "
              "- all frames kept (group statistics unreliable)")
        return
    for fr, reason, _kind in rejected[:max_listed]:
        print(f"  {'':<13}  REJECTED {fr['file'].name}: {reason}")
    if len(rejected) > max_listed:
        counts = Counter(kind for _fr, _reason, kind in rejected)
        breakdown = ", ".join(f"{k}: {v}" for k, v in sorted(counts.items()))
        print(f"  {'':<13}  ... and {len(rejected) - max_listed} more "
              f"({breakdown})")


def quality_control(frames):
    """
    Flag frames whose median level or noise is inconsistent with the group.
    Returns (accepted, rejected, info): rejected holds (frame, reason, kind)
    with kind in {"level", "noise"}, and info carries the thresholds that
    were applied (for reporting) plus the safety_net flag.
    """
    medians = np.array([f["median"] for f in frames])
    stds = np.array([f["std"] for f in frames])
    center = float(np.median(medians))
    spread = float(mad_std(medians))
    median_tol = MEDIAN_SIGMA * spread + MEDIAN_FLOOR_ADU
    std_center = float(np.median(stds))

    info = {
        "center": center,
        "spread": spread,
        "median_tol": median_tol,
        "std_center": std_center,
        "std_lo": std_center / STD_RATIO_MAX if std_center > 0 else None,
        "std_hi": std_center * STD_RATIO_MAX if std_center > 0 else None,
        "safety_net": False,
    }

    accepted, rejected = [], []
    for f in frames:
        reason = kind = None
        if abs(f["median"] - center) > median_tol:
            kind = "level"
            reason = (f"median {f['median']:.2f} deviates from group "
                      f"{center:.2f} by more than {median_tol:.2f} ADU")
        elif std_center > 0 and not (
            std_center / STD_RATIO_MAX <= f["std"] <= std_center * STD_RATIO_MAX
        ):
            kind = "noise"
            reason = (f"noise {f['std']:.2f} ADU is outside "
                      f"[{std_center / STD_RATIO_MAX:.2f}, "
                      f"{std_center * STD_RATIO_MAX:.2f}] (group: {std_center:.2f})")
        if reason:
            rejected.append((f, reason, kind))
        else:
            accepted.append(f)

    # Safety net: if QC would reject most of the group, the "group" statistics
    # themselves are unreliable - keep everything and just flag it.
    if len(accepted) < max(2, len(frames) // 2):
        info["safety_net"] = True
        return frames, [], info
    return accepted, rejected, info


def estimate_pair_noise(paths, max_pairs=MAX_RN_PAIRS):
    """
    Noise from differences of consecutive frame pairs:
    robust_std(a - b) / sqrt(2). For biases this is the read noise; for
    darks it also contains dark-current shot noise. Uses at most max_pairs
    pairs, loading two frames at a time. Returns the noise in ADU or None.
    """
    if len(paths) < 2:
        return None
    estimates = []
    prev = fits.getdata(paths[0], ext=0).astype(np.float32)
    for f in paths[1:1 + max_pairs]:
        cur = fits.getdata(f, ext=0).astype(np.float32)
        estimates.append(mad_std(prev - cur) / np.sqrt(2.0))
        prev = cur
    return float(np.median(estimates))


# ---------------------------------------------------------------------------
# Optional light-frame cleaning: hot-pixel repair and cosmic-ray rejection
# ---------------------------------------------------------------------------

def make_hot_pixel_mask(data, sigma=HOT_PIXEL_SIGMA):
    """
    Flag hot pixels from a master dark: pixels above
    clipped_median + sigma * clipped_std. Sigma-clipped statistics keep the
    threshold robust even when the pixel distribution is highly quantized
    (MAD near 0). Returns a bool array (True = hot pixel).
    """
    arr = np.asarray(data, dtype=np.float64)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return np.zeros(arr.shape, dtype=bool)
    _, med, std = sigma_clipped_stats(finite, sigma=3.0, maxiters=5)
    if std < 1e-6:
        std = float(np.std(finite))
    return arr > (float(med) + sigma * float(std))


def interpolate_bad_pixels(data, mask):
    """
    Replace each masked pixel with the median of its valid 3x3 neighbours.
    Pixels with no valid neighbour are left unchanged. Returns a float32 copy.

    The eight neighbours are gathered only at the masked positions - an
    (n_masked, 8) array - so the cost scales with the number of bad pixels
    instead of with a Python loop over them. Masked neighbours become NaN and
    drop out of the median, exactly as excluding them would; pixels outside
    the frame do the same via the NaN border. There is no order dependency
    (masked neighbours never contribute), so this matches the per-pixel loop
    it replaces bit for bit.
    """
    out = np.asarray(data, dtype=np.float32).copy()
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return out

    rows, cols = np.nonzero(mask)
    # NaN out the bad pixels, then pad by one so every 3x3 window is in range.
    vals = np.where(mask, np.nan, out).astype(np.float32)
    padded = np.pad(vals, 1, mode="constant", constant_values=np.nan)
    r, c = rows + 1, cols + 1

    neighbours = np.empty((rows.size, 8), dtype=np.float32)
    k = 0
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue                      # exclude the pixel itself
            neighbours[:, k] = padded[r + dr, c + dc]
            k += 1

    with warnings.catch_warnings():
        # all-NaN rows (a bad pixel with no valid neighbour) are expected
        warnings.simplefilter("ignore", RuntimeWarning)
        medians = np.nanmedian(neighbours, axis=1)

    usable = np.isfinite(medians)
    out[rows[usable], cols[usable]] = medians[usable]
    return out


def reject_cosmic_rays(ccd, gain, readnoise, sigclip=CR_SIGCLIP,
                       objlim=CR_OBJLIM):
    """
    Run L.A.Cosmic (ccdproc.cosmicray_lacosmic) on a calibrated frame. The
    cleaned values replace the cosmic rays in place; the CR mask is then
    dropped so compact()/write keep the repaired values instead of turning
    them into NaN. Data stays in ADU (gain_apply=False). Returns
    (ccd, n_cosmic_ray_pixels). Needs the 'astroscrappy' package.
    """
    cleaned = ccdp.cosmicray_lacosmic(
        ccd, sigclip=float(sigclip), objlim=float(objlim),
        gain=float(gain), readnoise=float(readnoise),
        gain_apply=False, verbose=False)
    n_cr = int(np.count_nonzero(cleaned.mask)) if cleaned.mask is not None else 0
    cleaned.mask = None
    return cleaned, n_cr
