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
from pathlib import Path

import numpy as np
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
    """
    try:
        return CCDData.read(path)
    except ValueError:
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

def scan_frames(paths):
    """
    Read frames one at a time and collect per-frame statistics.
    Returns a list of dicts: file, shape, median, std.
    """
    frames = []
    for f in paths:
        try:
            data = fits.getdata(f, ext=0).astype(np.float32)
        except Exception as exc:
            print(f"  WARNING: cannot read data of {Path(f).name}: {exc}")
            continue
        frames.append({
            "file": Path(f),
            "shape": data.shape,
            "median": float(np.median(data)),
            "std": float(mad_std(data)),
        })
    return frames


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


def print_frame_table(frames):
    """Print the per-frame statistics table (median / noise)."""
    name_w = max(len(fr["file"].name) for fr in frames)
    print(f"  {'frame':<{name_w}} {'median':>10} {'noise':>8}")
    for fr in frames:
        print(f"  {fr['file'].name:<{name_w}} {fr['median']:>10.2f} "
              f"{fr['std']:>8.2f}")


def quality_control(frames, skip_median_check=False):
    """
    Flag frames whose median level or noise is inconsistent with the group.
    With skip_median_check=True (level_align mode) only the noise check
    runs - level deviations are corrected by alignment instead of rejected.
    Returns (accepted, rejected) where rejected is a list of (frame, reason).
    """
    medians = np.array([f["median"] for f in frames])
    stds = np.array([f["std"] for f in frames])
    center = np.median(medians)
    spread = mad_std(medians)
    median_tol = MEDIAN_SIGMA * spread + MEDIAN_FLOOR_ADU
    std_center = np.median(stds)

    accepted, rejected = [], []
    for f in frames:
        reason = None
        if not skip_median_check and abs(f["median"] - center) > median_tol:
            reason = (f"median {f['median']:.2f} deviates from group "
                      f"{center:.2f} by more than {median_tol:.2f} ADU")
        elif std_center > 0 and not (
            std_center / STD_RATIO_MAX <= f["std"] <= std_center * STD_RATIO_MAX
        ):
            reason = (f"noise {f['std']:.2f} ADU is outside "
                      f"[{std_center / STD_RATIO_MAX:.2f}, "
                      f"{std_center * STD_RATIO_MAX:.2f}] (group: {std_center:.2f})")
        if reason:
            rejected.append((f, reason))
        else:
            accepted.append(f)

    # Safety net: if QC would reject most of the group, the "group" statistics
    # themselves are unreliable - keep everything and just warn.
    if len(accepted) < max(2, len(frames) // 2):
        print("  WARNING: QC would reject most frames - keeping all "
              "(group statistics unreliable)")
        return frames, []
    return accepted, rejected


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
    """
    out = np.asarray(data, dtype=np.float32).copy()
    if not mask.any():
        return out
    h, w = out.shape
    for r, c in zip(*np.where(mask)):
        r0, r1 = max(0, r - 1), min(h, r + 2)
        c0, c1 = max(0, c - 1), min(w, c + 2)
        patch = out[r0:r1, c0:c1]
        pmask = mask[r0:r1, c0:c1].copy()
        pmask[r - r0, c - c0] = True          # exclude the pixel itself
        valid = patch[~pmask]
        if valid.size:
            out[r, c] = float(np.median(valid))
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
