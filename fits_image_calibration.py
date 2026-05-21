#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CCD calibration utilities for building master calibration frames
(BIAS, DARK, FLAT) and calibrating science frames.

Features
--------
- Robust grouping of FITS frames by header keywords:
    (IMAGETYP, (XBINNING, YBINNING), GAIN, EXPTIME, rounded CCD-TEMP, FILTER).
- Temperature rounding with configurable step and tolerance reuse.
- Master selection with exposure- and temperature-matching tolerances.
- Sigma-clipped median combination with MAD-based dispersion.
- Optional use of existing masters for pre-correction.
- Science calibration pipeline: [bias] -> dark -> flat -> hot pixels -> cosmic rays.
- Hot pixel detection via sigma_clipped_stats; 3×3 neighbour median interpolation.
- Cosmic ray rejection via L.A.Cosmic (ccdproc.cosmicray_lacosmic).
- I/O in float32 to reduce memory footprint.
- All metadata read exclusively from FITS headers (no filename parsing).

Dependencies
------------
- numpy
- astropy (io.fits, units, stats)
- ccdproc (CCDData, combine, subtract_bias, subtract_dark, flat_correct,
           cosmicray_lacosmic)

This module is intended to be imported and called from a separate driver
script or notebook.

Author: Bogdan Alexandru Dumitru
License: MIT
"""

from __future__ import annotations

import logging
import math
from collections import defaultdict
from pathlib import Path
from typing import (
    Any,
    DefaultDict,
    Dict,
    Iterable,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.stats import mad_std
from ccdproc import (
    CCDData, combine, subtract_bias, subtract_dark, flat_correct,
    cosmicray_lacosmic,
)

logger = logging.getLogger(__name__)

# Key:
# (IMAGETYP, (XBINNING, YBINNING), GAIN, EXPTIME, CCD-TEMP-rounded, FILTER)
NumberKey = Tuple[
    str,
    Tuple[Optional[int], Optional[int]],
    Optional[int],
    Optional[float],
    Optional[float],
    Optional[str],   # FILTER — populated for FLAT frames; None for BIAS/DARK
]
GroupMap = Dict[NumberKey, List[Tuple[Path, fits.Header]]]

__all__ = [
    "load_files",
    "grouping",
    "create_master_file",
    "calibrate_files",
]


def _enable_debug_logger() -> None:
    """
    Attach a single StreamHandler with DEBUG level for this module logger.

    Designed for on-demand verbose diagnostics without affecting root logger.
    """
    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        logger.addHandler(handler)

    logger.setLevel(logging.DEBUG)
    logger.propagate = False


def load_files(
    directory: Union[str, Path],
    extensions: Sequence[str] = ("*.fit", "*.fits", "*.FIT", "*.FITS"),
) -> List[Path]:
    """
    Recursively collect FITS files from a directory.

    Parameters
    ----------
    directory : str or Path
        Root directory.
    extensions : sequence of str
        Filename patterns to match.

    Returns
    -------
    list of Path
        Sorted unique list of matching files.
    """
    root = Path(directory)
    if not root.is_dir():
        return []

    out: List[Path] = []
    for pattern in extensions:
        out.extend(root.rglob(pattern))

    # Use set to avoid duplicates from overlapping patterns.
    return sorted(set(out))


def _safe_int(value: Any) -> Optional[int]:
    """Convert value to nearest int, or return None on failure."""
    try:
        return int(round(float(value))) if value is not None else None
    except (TypeError, ValueError):
        return None


def _safe_float(value: Any) -> Optional[float]:
    """Convert value to float, or return None on failure."""
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def _round_temp(value: Any, step: float = 0.5) -> Optional[float]:
    """
    Round CCD temperature to a grid defined by 'step'.

    Returns None if conversion fails.
    """
    v = _safe_float(value)
    if v is None:
        return None
    return round(v / step) * step


def _format_bin(binning: Tuple[Optional[int], Optional[int]]) -> str:
    """Format binning tuple for filenames."""
    x, y = binning
    return f"b{('NA' if x is None else x)}x{('NA' if y is None else y)}"


def _float_str(v: Optional[float], fmt: str = ".3g") -> str:
    """Format optional float for filenames."""
    return f"{v:{fmt}}" if v is not None else "NA"


def _get_object_name(file_path: Path, base_dir: Path) -> Optional[str]:
    """
    Return the immediate subdirectory name of file_path relative to base_dir.

    E.g.  raw/LIGHT/NGC1234/img.fits  ->  'NGC1234'
    Returns None if the file is directly inside base_dir.
    """
    try:
        rel = file_path.relative_to(base_dir)
        if len(rel.parts) > 1:
            return rel.parts[0]
    except ValueError:
        pass
    return None


def grouping(
    file_list: Iterable[Union[Path, str]],
    temp_tolerance: float = 2.0,
    round_step: float = 0.5,
) -> GroupMap:
    """
    Group FITS files based on header keywords.

    Group key:
        (IMAGETYP, (XBINNING, YBINNING), GAIN, EXPTIME, rounded CCD-TEMP, FILTER)

    FILTER handling:
        - All metadata is read exclusively from FITS headers; filenames are
          never parsed and files are never modified by this function.
        - For FLAT/LIGHT frames: FILTER keyword is used as-is; defaults to
          "Clear" when the keyword is absent.
        - For BIAS and DARK, FILTER is always None in the key.

    CCD-TEMP handling:
        - Temperature is rounded using 'round_step'.
        - If an existing group has compatible metadata and
            |ΔT| <= temp_tolerance, the file is attached to that group.

    Only primary headers are read for performance.

    Parameters
    ----------
    file_list : iterable of Path or str
        Input FITS file paths.
    temp_tolerance : float
        Absolute temperature tolerance [deg C] for reusing existing groups.
    round_step : float
        Rounding step for CCD-TEMP.

    Returns
    -------
    dict
        Mapping from grouping key to list of (Path, Header).
    """
    files = [Path(f) for f in file_list]
    groups: DefaultDict[NumberKey, List[Tuple[Path, fits.Header]]] = defaultdict(list)

    for file in files:
        if not file.exists():
            logger.warning("File not found: %s", file)
            continue

        try:
            header = fits.getheader(file, ext=0)

            filetype = str(header.get("IMAGETYP", "")).strip().upper()
            binning = (
                _safe_int(header.get("XBINNING")),
                _safe_int(header.get("YBINNING")),
            )
            gain = _safe_int(header.get("GAIN"))
            exposure = _safe_float(header.get("EXPTIME"))
            temperature = _round_temp(header.get("CCD-TEMP"), step=round_step)

            # FILTER is read exclusively from the FITS header.
            # BIAS and DARK frames never carry a filter value.
            # For FLAT/LIGHT frames, default to "Clear" when the keyword is absent.
            filter_name: Optional[str] = None
            if filetype not in {"BIAS", "DARK"}:
                header_filter = str(header.get("FILTER", "")).strip()
                filter_name = header_filter if header_filter else "Clear"

            matched_key: Optional[NumberKey] = None

            # Search for an existing compatible key within temperature tolerance.
            for k_type, k_bin, k_gain, k_exp, k_temp, k_filter in list(groups.keys()):
                if k_type != filetype:
                    continue
                if k_bin != binning:
                    continue
                if k_gain != gain:
                    continue
                if k_filter != filter_name:
                    continue

                # Exposure must match (if both present) for DARK; for others we
                # only require exact equality when defined.
                if (k_exp is None) != (exposure is None):
                    continue
                if k_exp is not None and exposure is not None:
                    if not math.isclose(k_exp, exposure, rel_tol=1e-6):
                        continue

                if k_temp is None or temperature is None:
                    continue

                if abs(k_temp - temperature) <= temp_tolerance:
                    matched_key = (k_type, k_bin, k_gain, k_exp, k_temp, k_filter)
                    break

            key: NumberKey = matched_key or (
                filetype,
                binning,
                gain,
                exposure,
                temperature,
                filter_name,
            )
            groups[key].append((file, header))

        except Exception:
            logger.exception("Failed to read header for %s", file)
            continue

    logger.info("Grouped %d files into %d groups", len(files), len(groups))
    return dict(groups)


def _find_single_match(
    master_groups: GroupMap,
    want_type: str,
    binning: Tuple[Optional[int], Optional[int]],
    gain: Optional[int],
    exposure: Optional[float],
    temp: Optional[float],
    filter_name: Optional[str] = None,
    tol_exp_rel: float = 1e-6,
    tol_temp_abs: float = 2.0,
) -> Optional[Tuple[Path, fits.Header]]:
    """
    Select a single compatible master frame.

    Matching rules
    --------------
    - IMAGETYP == want_type
    - Exact binning and gain
    - For DARK: exposure must match within tol_exp_rel
    - For BIAS/FLAT: exposure not constrained (unique solution assumed)
    - For FLAT: filter_name must match exactly
    - If both master and target have CCD-TEMP: |ΔT| <= tol_temp_abs

    If multiple keys are compatible, raise RuntimeError to avoid ambiguity.
    If multiple files share the same key, choose the newest by mtime.

    Returns
    -------
    (Path, Header) or None
    """
    candidate_keys: List[NumberKey] = []

    for m_type, m_bin, m_gain, m_exp, m_temp, m_filter in master_groups.keys():
        if m_type != want_type:
            continue
        if m_bin != binning:
            continue
        if m_gain != gain:
            continue

        if want_type == "DARK":
            if m_exp is None or exposure is None:
                continue
            if not math.isclose(m_exp, exposure, rel_tol=tol_exp_rel):
                continue

        if want_type == "FLAT":
            if m_filter != filter_name:
                continue

        if m_temp is not None and temp is not None:
            if abs(m_temp - temp) > tol_temp_abs:
                continue

        candidate_keys.append((m_type, m_bin, m_gain, m_exp, m_temp, m_filter))

    if not candidate_keys:
        return None

    if len(candidate_keys) > 1:
        raise RuntimeError(
            f"Multiple master {want_type} keys for "
            f"bin={binning} gain={gain} exp={exposure} temp={temp} filter={filter_name}"
        )

    key = candidate_keys[0]
    files = master_groups[key]

    if len(files) == 1:
        return files[0]

    # Prefer the newest master file for this key.
    return sorted(files, key=lambda t: t[0].stat().st_mtime, reverse=True)[0]


def _find_single_match_interactive(
    master_groups: GroupMap,
    want_type: str,
    binning: Tuple[Optional[int], Optional[int]],
    gain: Optional[int],
    exposure: Optional[float],
    temp: Optional[float],
    filter_name: Optional[str] = None,
    tol_exp_rel: float = 1e-6,
    tol_temp_abs: float = 2.0,
) -> Optional[Tuple[Path, fits.Header]]:
    """
    Wrapper over _find_single_match that handles ambiguous matches interactively.

    If multiple compatible masters exist, lists them and asks the user
    to pick one via stdin. Blocks until a valid choice is entered.

    Returns
    -------
    (Path, Header) or None
    """
    try:
        return _find_single_match(
            master_groups, want_type, binning, gain, exposure, temp,
            filter_name, tol_exp_rel, tol_temp_abs,
        )
    except RuntimeError:
        # Collect all candidate keys manually to show the user
        candidate_files: List[Tuple[Path, fits.Header]] = []
        for key, file_list in master_groups.items():
            m_type, m_bin, m_gain, m_exp, m_temp, m_filter = key
            if m_type != want_type:
                continue
            if m_bin != binning:
                continue
            if m_gain != gain:
                continue
            if want_type == "DARK":
                if m_exp is None or exposure is None:
                    continue
                if not math.isclose(m_exp, exposure, rel_tol=tol_exp_rel):
                    continue
            if want_type == "FLAT":
                if m_filter != filter_name:
                    continue
            if m_temp is not None and temp is not None:
                if abs(m_temp - temp) > tol_temp_abs:
                    continue
            candidate_files.extend(file_list)

        print(
            f"\n[AMBIGUOUS] Multiple master {want_type} files found "
            f"for bin={binning} gain={gain} exp={exposure} temp={temp} filter={filter_name}:"
        )
        for i, (fpath, _) in enumerate(candidate_files):
            print(f"  [{i + 1}] {fpath}")

        while True:
            try:
                raw = input(f"Select master {want_type} [1-{len(candidate_files)}]: ").strip()
                idx = int(raw) - 1
                if 0 <= idx < len(candidate_files):
                    chosen = candidate_files[idx]
                    print(f"  -> Using: {chosen[0]}")
                    return chosen
                else:
                    print(f"  Invalid choice. Enter a number between 1 and {len(candidate_files)}.")
            except ValueError:
                print(f"  Invalid input. Enter a number between 1 and {len(candidate_files)}.")
            except EOFError:
                chosen = candidate_files[0]
                logger.warning(
                    "EOF received during master selection, falling back to: %s", chosen[0]
                )
                return chosen


def _normalize_master_paths(
    master_files: Optional[Union[str, Path, Iterable[Union[str, Path]]]],
) -> List[Path]:
    """
    Normalize user-provided master paths.

    - If None: return empty list.
    - If file: return [file].
    - If directory: recursively collect FITS files.
    - If iterable: convert each element to Path.
    """
    if master_files is None:
        return []

    # Single path-like
    if isinstance(master_files, (str, Path)):
        p = Path(master_files)
        if p.is_dir():
            return load_files(p, ("*.fit", "*.fits"))
        if p.is_file() and p.suffix.lower() in (".fit", ".fits"):
            return [p]
        return []

    # Iterable of paths: expand dirs, keep only FITS files
    out: List[Path] = []
    for item in master_files:
        p = Path(item)
        if p.is_dir():
            out.extend(load_files(p, ("*.fit", "*.fits")))
        elif p.is_file() and p.suffix.lower() in (".fit", ".fits"):
            out.append(p)
    return out


def _add_history(meta: Any, text: str) -> None:
    """
    Append a HISTORY entry to FITS-like metadata.

    Supports both astropy.io.fits.Header (via add_history)
    and dict-like containers.
    """
    if hasattr(meta, "add_history"):
        meta.add_history(text)
        return

    prev = meta.get("HISTORY")
    if prev is None:
        meta["HISTORY"] = [text]
    elif isinstance(prev, list):
        prev.append(text)
    else:
        meta["HISTORY"] = [str(prev), text]


def _inv_median(ccd_obj: CCDData) -> float:
    """Scale factor for flat combination: 1 / median(frame)."""
    arr = np.ma.array(ccd_obj.data)
    med = np.ma.median(arr)
    if np.ma.is_masked(med) or float(med) == 0.0:
        return 1.0
    return float(1.0 / med)


def _get_fits_unit(header: Any) -> u.Unit:
    """
    Determine the physical unit from a FITS header BUNIT keyword.

    Falls back to ADU if BUNIT is absent, empty, or unrecognized by astropy.

    Parameters
    ----------
    header : fits.Header or dict-like
        Primary HDU header.

    Returns
    -------
    astropy.units.Unit
    """
    bunit = None
    if hasattr(header, "get"):
        bunit = header.get("BUNIT")

    if bunit:
        bunit_str = str(bunit).strip()
        if bunit_str:
            try:
                return u.Unit(bunit_str)
            except ValueError:
                logger.warning(
                    "Unrecognized BUNIT '%s', falling back to ADU", bunit_str
                )

    return u.adu


def _to_ccd(hdu: Union[fits.HDUList, fits.hdu.image.PrimaryHDU]) -> CCDData:
    """
    Convert a PrimaryHDU or HDUList to CCDData in float32.

    The physical unit is read from the BUNIT header keyword.
    Falls back to ADU if BUNIT is absent or unrecognized.

    If HDUList is provided, the primary HDU is used.
    """
    primary = hdu[0] if isinstance(hdu, fits.HDUList) else hdu
    data = getattr(primary, "data", None)
    if data is None:
        raise ValueError("HDU data is None")
    header = getattr(primary, "header", {})
    unit = _get_fits_unit(header)
    return CCDData(
        data.astype(np.float32, copy=False),
        meta=header,
        unit=unit,
    )


def _make_hot_pixel_mask(
    master_dark: CCDData,
    sigma: float = 5.0,
) -> np.ndarray:
    """
    Derive a bad-pixel mask from a master dark frame.

    Uses sigma-clipped statistics so the threshold is robust even when the
    pixel distribution is highly quantized (MAD = 0).

    Pixels whose value exceeds  clipped_median + sigma * clipped_std  are
    flagged as hot.

    Returns
    -------
    numpy bool array  (True = hot pixel)
    """
    from astropy.stats import sigma_clipped_stats

    data = np.asarray(master_dark.data, dtype=np.float64)
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return np.zeros(data.shape, dtype=bool)

    _, med_clipped, std_clipped = sigma_clipped_stats(finite, sigma=3.0, maxiters=5)
    # Guard against degenerate case where clipped std is still near zero
    if std_clipped < 1e-6:
        std_clipped = float(np.std(finite))
    threshold = float(med_clipped) + sigma * float(std_clipped)
    mask      = data > threshold
    logger.debug(
        "Hot pixel mask: median=%.3f  clipped_std=%.4f  threshold=%.3f  "
        "flagged=%d / %d (%.3f%%)",
        med_clipped, std_clipped, threshold,
        int(mask.sum()), mask.size,
        100.0 * mask.sum() / mask.size,
    )
    return mask


def _interpolate_bad_pixels(
    data: np.ndarray,
    mask: np.ndarray,
) -> np.ndarray:
    """
    Replace each masked pixel with the median of its valid 3×3 neighbours.

    Pixels with no valid neighbour are left unchanged.

    Returns
    -------
    float32 array with bad pixels interpolated in-place copy.
    """
    if not mask.any():
        return data.astype(np.float32, copy=True)

    out = data.astype(np.float32, copy=True)
    h, w = data.shape
    bad_rows, bad_cols = np.where(mask)

    for r, c in zip(bad_rows.tolist(), bad_cols.tolist()):
        r0, r1 = max(0, r - 1), min(h, r + 2)
        c0, c1 = max(0, c - 1), min(w, c + 2)
        patch      = out[r0:r1, c0:c1]
        patch_mask = mask[r0:r1, c0:c1].copy()
        patch_mask[r - r0, c - c0] = True          # exclude the pixel itself
        valid = patch[~patch_mask]
        if valid.size > 0:
            out[r, c] = float(np.median(valid))

    return out


def _apply_cosmicray_rejection(
    ccd:       CCDData,
    gain:      float = 1.0,
    readnoise: float = 8.0,
    sigclip:   float = 4.5,
    objlim:    float = 5.0,
) -> Tuple[CCDData, int]:
    """
    Run L.A.Cosmic cosmic-ray detection on a calibrated frame.

    Parameters
    ----------
    ccd       : calibrated CCDData (bias+dark+flat corrected)
    gain      : CCD gain [e⁻/ADU]  — read from EGAIN header keyword (falls back to GAIN)
    readnoise : read noise [e⁻]    — read from RDNOISE or config default
    sigclip   : detection threshold in sigma above background
    objlim    : minimum contrast between CR and underlying source

    Returns
    -------
    (cleaned_ccd, n_cr_pixels)
        cleaned_ccd  — CCDData with CR pixels replaced and mask set
        n_cr_pixels  — number of pixels flagged as cosmic rays
    """
    cleaned = cosmicray_lacosmic(
        ccd,
        sigclip=float(sigclip),
        objlim=float(objlim),
        gain=float(gain),
        readnoise=float(readnoise),
        gain_apply=False,   # keep data in ADU
        verbose=False,
    )
    cr_mask = cleaned.mask
    n_cr    = int(np.sum(cr_mask)) if cr_mask is not None else 0
    return cleaned, n_cr


def create_master_file(
    dir_path: Union[str, Path],
    file_type: str,
    output: Union[str, Path],
    master_files: Optional[Union[str, Path, Iterable[Union[str, Path]]]] = None,
    no_bias: bool = False,
    dark_scale: bool = False,
    debug: bool = False,
    mem_limit: Optional[int] = None,
) -> List[Path]:
    """
    Build master calibration frames for data in ``dir_path/<file_type>``.

    For each compatible group:
    - BIAS, DARK:
        median combine with sigma clipping.
    - FLAT:
        median combine with per-frame scaling by 1 / median(level),
        optional pre-subtraction of master bias and master dark.

    Existing masters
    ----------------
    - If ``master_files`` is provided, they are used as candidates
        for pre-correction (bias/dark) of input frames.
    - If not provided and file_type is DARK or FLAT, existing masters
        are searched in the ``output`` directory.

    Parameters
    ----------
    dir_path : str or Path
        Root directory containing subdirectories named by IMAGETYP.
    file_type : str
        One of: 'BIAS', 'DARK', 'FLAT'.
    output : str or Path
        Directory where master files will be written.
    master_files : path-like or iterable, optional
        Paths to existing master frames or directory containing them.
    no_bias : bool
        If True, skip bias subtraction when building DARK/FLAT masters.
    dark_scale : bool
        If True, scale darks by exposure time in subtract_dark.
    debug : bool
        If True, enable verbose logging for this module.
    mem_limit : int, optional
        Memory limit in bytes for ccdproc.combine (if supported).

    Returns
    -------
    list of Path
        Paths to successfully created master files.
    """
    if debug:
        _enable_debug_logger()

    dir_path = Path(dir_path)
    output = Path(output)
    file_type_norm = str(file_type).strip().upper()

    # Collect input frames from subdirectory named after file_type.
    files = load_files(dir_path / file_type_norm, ("*.fit", "*.fits"))
    groups = grouping(files)

    # Prepare master candidates (for pre-correction).
    master_paths = _normalize_master_paths(master_files)
    if not master_paths and file_type_norm in {"DARK", "FLAT"}:
        master_paths = load_files(output, ("*.fit", "*.fits"))
    master_groups = grouping(master_paths) if master_paths else {}

    created: List[Path] = []

    for (ftype, binning, gain, exposure, temperature, filter_name), file_list in groups.items():
        # Only process frames whose IMAGETYP matches requested file_type.
        if ftype != file_type_norm:
            continue
        if not file_list:
            continue

        master_bias: Optional[CCDData] = None
        master_dark: Optional[CCDData] = None

        # Try to load matching masters for pre-correction of input frames.
        if master_groups:
            try:
                if file_type_norm in {"DARK", "FLAT"} and not no_bias:
                    bias_match = _find_single_match(
                        master_groups, "BIAS", binning, gain, None, None
                    )
                    if bias_match:
                        bias_path, _ = bias_match
                        with fits.open(bias_path) as hdul:
                            master_bias = _to_ccd(hdul)

                if file_type_norm == "FLAT":
                    dark_match = _find_single_match(
                        master_groups, "DARK", binning, gain, exposure, temperature
                    )
                    if dark_match:
                        dark_path, _ = dark_match
                        with fits.open(dark_path) as hdul:
                            master_dark = _to_ccd(hdul)
            except RuntimeError:
                logger.exception(
                    "Master file matching error for group %s", (ftype, binning, gain)
                )

        ccd_list: List[CCDData] = []

        # Prepare CCDData objects with optional bias/dark correction.
        for fpath, _ in file_list:
            try:
                with fits.open(fpath) as hdul:
                    ccd = _to_ccd(hdul)

                if file_type_norm in {"DARK", "FLAT"}:
                    if master_bias is not None:
                        ccd = subtract_bias(ccd, master_bias)
                        _add_history(ccd.meta, "Master bias subtracted")
                    elif not no_bias:
                        logger.warning("Bias subtraction skipped for %s", fpath.name)

                if file_type_norm == "FLAT":
                    if master_dark is not None:
                        ccd = subtract_dark(
                            ccd,
                            master_dark,
                            exposure_time="EXPTIME",
                            exposure_unit=u.second,
                            scale=dark_scale,
                        )
                        _add_history(ccd.meta, "Master dark subtracted")
                    else:
                        logger.warning("Dark subtraction skipped for %s", fpath.name)

                ccd_list.append(ccd)

            except Exception:
                logger.exception("Failed to prepare CCDData from %s", fpath)

        if not ccd_list:
            logger.info("No valid frames to combine for group %s %s", ftype, binning)
            continue

        master: Optional[CCDData] = None

        # Configure combine() keyword arguments.
        combine_kwargs: Dict[str, Any] = dict(
            method="median",
            sigma_clip=True,
            sigma_clip_low_thresh=3,
            sigma_clip_high_thresh=3,
            sigma_clip_func=np.ma.median,
            sigma_clip_dev_func=mad_std,
        )
        if mem_limit is not None:
            combine_kwargs["mem_limit"] = float(mem_limit)

        if file_type_norm in {"BIAS", "DARK"}:
            master = combine(ccd_list, **combine_kwargs)
            if master is not None:
                _add_history(
                    master.meta,
                    f"Master {file_type_norm} bin={binning} "
                    f"gain={gain} exp={_float_str(exposure)}s temp={temperature}C",
                )
                _add_history(master.meta, f"Combined from {len(ccd_list)} files")

        elif file_type_norm == "FLAT":
            combine_kwargs_flat: Dict[str, Any] = dict(
                method="median",
                scale=_inv_median,
                sigma_clip=True,
                sigma_clip_low_thresh=3,
                sigma_clip_high_thresh=3,
                sigma_clip_func=np.ma.median,
                sigma_clip_dev_func=mad_std,
            )
            if mem_limit is not None:
                combine_kwargs_flat["mem_limit"] = float(mem_limit)

            master = combine(ccd_list, **combine_kwargs_flat)
            if master is not None:
                _add_history(
                    master.meta,
                    f"Master FLAT bin={binning} gain={gain} "
                    f"exp={_float_str(exposure)}s temp={temperature}C "
                    f"filter={filter_name}",
                )
                _add_history(master.meta, f"Combined from {len(ccd_list)} files")

        if master is None:
            continue

        output.mkdir(parents=True, exist_ok=True)

        bin_str = _format_bin(binning)
        gain_str = f"g{gain}" if gain is not None else "gNA"
        exp_str = f"{_float_str(exposure)}s" if exposure is not None else "sNA"
        temp_str = f"{_float_str(temperature)}C" if temperature is not None else "CNA"
        filter_str = f"_{filter_name}" if filter_name is not None else ""

        out_name = (
            output
            / f"master_{file_type_norm.lower()}_{bin_str}_{gain_str}_{exp_str}_{temp_str}{filter_str}.fits"
        )

        try:
            master.data = master.data.astype(np.float32, copy=False)
            primary_hdu = fits.PrimaryHDU(
                data=master.data,
                header=fits.Header(master.meta),
            )
            fits.HDUList([primary_hdu]).writeto(out_name, overwrite=True)
            logger.info("Created master file %s from %d files", out_name, len(ccd_list))
            created.append(out_name)
        except Exception:
            logger.exception("Failed writing master to %s", out_name)
    return created


def calibrate_files(
    input_dir: Union[str, Path],
    file_type: str,
    output_dir: Union[str, Path],
    master_files: Union[str, Path, Iterable[Union[str, Path]]],
    no_bias: bool = False,
    dark_scale: bool = False,
    debug: bool = False,
    fix_hot_pixels:   bool  = True,
    hot_pixel_sigma:  float = 5.0,
    fix_cosmic_rays:  bool  = True,
    cr_sigclip:       float = 4.5,
    cr_objlim:        float = 5.0,
    cr_readnoise:     float = 8.0,
) -> List[Path]:
    """
    Calibrate frames from ``input_dir/<file_type>/`` using master frames.

    Science frames are discovered recursively under ``input_dir/<file_type>/``.
    The immediate subdirectory name is used as the object name.

    Output layout
    -------------
    ``output_dir/<object>/<filter>/filename_cal.fits``

    Calibration sequence
    --------------------
    For each group:
        1. Subtract master bias (if available and not no_bias).
        2. Subtract master dark (matched in exposure and temperature).
        3. Apply flat-field correction.
        4. Interpolate hot pixels detected from master dark (if fix_hot_pixels).
        5. Reject cosmic rays using L.A.Cosmic (if fix_cosmic_rays).

    Parameters
    ----------
    input_dir : str or Path
        Root directory containing ``<file_type>/<object>/`` subdirectories.
    file_type : str
        Science image type to calibrate (e.g. 'LIGHT').
    output_dir : str or Path
        Root output directory for calibrated FITS files.
    master_files : path-like or iterable
        Master calibration frames (or directory containing them).
    no_bias : bool
        If True, skip bias subtraction.
    dark_scale : bool
        If True, scale master darks by exposure time in subtract_dark.
    debug : bool
        If True, enable verbose logging.
    fix_hot_pixels : bool
        Detect hot pixels from master dark and interpolate them in each frame.
    hot_pixel_sigma : float
        Threshold multiplier for hot pixel detection: clipped_median + sigma * clipped_std.
    fix_cosmic_rays : bool
        Run L.A.Cosmic cosmic-ray rejection on each calibrated frame.
    cr_sigclip : float
        L.A.Cosmic detection threshold [sigma above background].
    cr_objlim : float
        L.A.Cosmic minimum contrast between CR and underlying source.
    cr_readnoise : float
        Default read noise [e⁻] used when RDNOISE is absent from the header.

    Returns
    -------
    list of Path
        Paths of calibrated files.
    """
    if debug:
        _enable_debug_logger()

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    file_type_norm = str(file_type).strip().upper()
    light_base = input_dir / file_type_norm

    files = load_files(light_base, ("*.fit", "*.fits"))
    groups = grouping(files)

    master_paths = _normalize_master_paths(master_files)
    master_groups = grouping(master_paths)

    created: List[Path] = []

    for (imagetyp, binning, gain, exposure, temperature, filter_name), file_list in groups.items():
        # Only treat frames whose header IMAGETYP matches requested file_type.
        if imagetyp != file_type_norm:
            continue
        if not file_list:
            continue

        filter_label = filter_name if filter_name is not None else "nofilter"

        # Find best master matches for this group.
        bias_match = (
            None
            if no_bias
            else _find_single_match_interactive(master_groups, "BIAS", binning, gain, None, None)
        )
        dark_match = _find_single_match_interactive(
            master_groups, "DARK", binning, gain, exposure, temperature
        )
        flat_match = _find_single_match_interactive(
            master_groups, "FLAT", binning, gain, None, None,
            filter_name=filter_name,
        )

        master_bias = None
        master_dark = None
        master_flat = None
        dark_path: Optional[Path] = None

        if bias_match:
            bias_path, _ = bias_match
            with fits.open(bias_path) as hdul:
                master_bias = _to_ccd(hdul)

        if dark_match:
            dark_path, _ = dark_match
            with fits.open(dark_path) as hdul:
                master_dark = _to_ccd(hdul)

        if flat_match:
            flat_path, _ = flat_match
            with fits.open(flat_path) as hdul:
                master_flat = _to_ccd(hdul)

        # Build hot pixel mask once per group from the matching master dark.
        hot_mask: Optional[np.ndarray] = None
        if fix_hot_pixels and master_dark is not None:
            hot_mask = _make_hot_pixel_mask(master_dark, sigma=hot_pixel_sigma)
            logger.info(
                "  Hot pixel mask: %d pixels (sigma=%.1f) from %s",
                int(hot_mask.sum()), hot_pixel_sigma,
                dark_path.name if dark_path is not None else "unknown",
            )

        for fp, _ in file_list:
            try:
                with fits.open(fp) as hdul:
                    ccd = _to_ccd(hdul)

                if master_bias is not None:
                    ccd = subtract_bias(ccd, master_bias)
                    _add_history(ccd.meta, "Master bias subtracted")
                elif not no_bias:
                    logger.warning("Missing BIAS for %s", fp)

                if master_dark is not None:
                    ccd = subtract_dark(
                        ccd,
                        master_dark,
                        exposure_time="EXPTIME",
                        exposure_unit=u.second,
                        scale=dark_scale,
                    )
                    _add_history(ccd.meta, "Master dark subtracted")
                else:
                    logger.warning("Missing DARK for %s", fp)

                if master_flat is not None:
                    flat_data = np.asarray(master_flat.data, dtype=float)
                    finite = flat_data[np.isfinite(flat_data)]
                    flat_floor = max(1e-6, 0.01 * float(np.median(finite))) if finite.size > 0 else 1e-6
                    ccd = flat_correct(ccd, master_flat, min_value=flat_floor, norm_value=1.0)
                    _add_history(ccd.meta, f"Flat-field corrected (min_value={flat_floor:.3g})")
                else:
                    logger.warning("Missing FLAT for %s", fp)

                # --- Hot pixel interpolation ---
                if hot_mask is not None and hot_mask.any():
                    ccd.data = _interpolate_bad_pixels(ccd.data, hot_mask)
                    _add_history(
                        ccd.meta,
                        f"Hot pixel correction: {int(hot_mask.sum())} pixels interpolated",
                    )

                # --- Cosmic ray rejection ---
                if fix_cosmic_rays:
                    try:
                        gain_val = max(0.1,
                                         _safe_float(ccd.meta.get("EGAIN"))
                                         or _safe_float(ccd.meta.get("GAIN"))
                                         or 1.0)
                        rn_val   = max(0.1, _safe_float(ccd.meta.get("RDNOISE")) or cr_readnoise)
                        ccd, n_cr = _apply_cosmicray_rejection(
                            ccd,
                            gain=gain_val,
                            readnoise=rn_val,
                            sigclip=cr_sigclip,
                            objlim=cr_objlim,
                        )
                        if n_cr > 0:
                            logger.info("  CR rejection: %d pixel(s) in %s", n_cr, fp.name)
                            _add_history(
                                ccd.meta,
                                f"Cosmic ray rejection: {n_cr} pixel(s) cleaned (L.A.Cosmic)",
                            )
                    except Exception:
                        logger.warning("  Cosmic ray rejection failed for %s", fp.name)

                # Build output path: calibrated/<object>/<filter>/stem_cal.fits
                obj_name = _get_object_name(fp, light_base)
                obj_label = obj_name if obj_name else "unknown"
                out_dir = output_dir / obj_label / filter_label
                out_dir.mkdir(parents=True, exist_ok=True)
                out_path = out_dir / f"{fp.stem}_cal.fits"

                ccd.data = ccd.data.astype(np.float32, copy=False)
                ccd.write(out_path, overwrite=True)
                created.append(out_path)

            except Exception:
                logger.exception("Calibration failed for %s", fp)

    return created
