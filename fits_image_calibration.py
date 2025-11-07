#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CCD calibration utilities for building master calibration frames
(BIAS, DARK, FLAT) and calibrating science frames.

Features
--------
- Robust grouping of FITS frames by:
  (IMAGETYP, (XBINNING, YBINNING), GAIN, EXPTIME, rounded CCD-TEMP).
- Temperature rounding with configurable step and tolerance reuse.
- Master selection with exposure- and temperature-matching tolerances.
- Sigma-clipped median combination with MAD-based dispersion.
- Optional use of existing masters for pre-correction.
- Science calibration pipeline: [bias] -> dark -> flat.
- I/O in float32 to reduce memory footprint.

Dependencies
------------
- numpy
- astropy (io.fits, units, stats)
- ccdproc (CCDData, combine, subtract_bias, subtract_dark, flat_correct)

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
from astropy.stats import SigmaClip, mad_std
from ccdproc import CCDData, combine, subtract_bias, subtract_dark, flat_correct

logger = logging.getLogger(__name__)

# Key:
# (IMAGETYP, (XBINNING, YBINNING), GAIN, EXPTIME, CCD-TEMP-rounded)
NumberKey = Tuple[str, Tuple[Optional[int], Optional[int]], Optional[int], Optional[float], Optional[float]]
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


def load_files(directory: Union[str, Path], extensions: Sequence[str] = ("*.fit", "*.fits")) -> List[Path]:
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


def grouping(
    file_list: Iterable[Union[Path, str]],
    temp_tolerance: float = 2.0,
    round_step: float = 0.5,
) -> GroupMap:
    """
    Group FITS files based on header keywords.

    Group key:
        (IMAGETYP, (XBINNING, YBINNING), GAIN, EXPTIME, rounded CCD-TEMP)

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
            binning = (_safe_int(header.get("XBINNING")), _safe_int(header.get("YBINNING")))
            gain = _safe_int(header.get("GAIN"))
            exposure = _safe_float(header.get("EXPTIME"))
            temperature = _round_temp(header.get("CCD-TEMP"), step=round_step)

            matched_key: Optional[NumberKey] = None

            # Search for an existing compatible key within temperature tolerance.
            for (k_type, k_bin, k_gain, k_exp, k_temp) in list(groups.keys()):
                if k_type != filetype:
                    continue
                if k_bin != binning:
                    continue
                if k_gain != gain:
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
                    matched_key = (k_type, k_bin, k_gain, k_exp, k_temp)
                    break

            key: NumberKey = matched_key or (filetype, binning, gain, exposure, temperature)
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
    - If both master and target have CCD-TEMP: |ΔT| <= tol_temp_abs

    If multiple keys are compatible, raise RuntimeError to avoid ambiguity.
    If multiple files share the same key, choose the newest by mtime.

    Returns
    -------
    (Path, Header) or None
    """
    candidate_keys: List[NumberKey] = []

    for (m_type, m_bin, m_gain, m_exp, m_temp) in master_groups.keys():
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

        if m_temp is not None and temp is not None:
            if abs(m_temp - temp) > tol_temp_abs:
                continue

        candidate_keys.append((m_type, m_bin, m_gain, m_exp, m_temp))

    if not candidate_keys:
        return None

    if len(candidate_keys) > 1:
        raise RuntimeError(
            f"Multiple master {want_type} keys for bin={binning} gain={gain}"
        )

    key = candidate_keys[0]
    files = master_groups[key]

    if len(files) == 1:
        return files[0]

    # Prefer the newest master file for this key.
    return sorted(files, key=lambda t: t[0].stat().st_mtime, reverse=True)[0]

def _normalize_master_paths(master_files: Optional[Union[str, Path, Iterable[Union[str, Path]]]]) -> List[Path]:
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


def _to_ccd(hdu: fits.hdu.image.PrimaryHDU) -> CCDData:
    """
    Convert a PrimaryHDU to CCDData in ADU with float32 data.

    Assumes the primary HDU contains the image data.
    """
    return CCDData(
        hdu.data.astype(np.float32, copy=False),
        meta=hdu.header,
        unit=u.adu,
    )


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

    sigma = SigmaClip(sigma=3.0, maxiters=5)

    for (ftype, binning, gain, exposure, temperature), file_list in groups.items():
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
                            master_bias = _to_ccd(hdul[0])

                if file_type_norm == "FLAT":
                    dark_match = _find_single_match(
                        master_groups, "DARK", binning, gain, exposure, temperature
                    )
                    if dark_match:
                        dark_path, _ = dark_match
                        with fits.open(dark_path) as hdul:
                            master_dark = _to_ccd(hdul[0])
            except RuntimeError:
                logger.exception(
                    "Master file matching error for group %s", (ftype, binning, gain)
                )

        ccd_list: List[CCDData] = []

        # Prepare CCDData objects with optional bias/dark correction.
        for fpath, _ in file_list:
            try:
                with fits.open(fpath) as hdul:
                    ccd = _to_ccd(hdul[0])

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
            sigma_clip=sigma,
            sigma_clip_func=np.ma.median,
            sigma_clip_dev_func=mad_std,
        )
        if mem_limit is not None:
            combine_kwargs["mem_limit"] = float(mem_limit)

        if file_type_norm in {"BIAS", "DARK"}:
            master = combine(ccd_list, **combine_kwargs)
            _add_history(
                master.meta,
                f"Master {file_type_norm} bin={binning} "
                f"gain={gain} exp={_float_str(exposure)}s temp={temperature}C",
            )
            _add_history(master.meta, f"Combined from {len(ccd_list)} files")

        elif file_type_norm == "FLAT":
            # Scale each flat frame by inverse median to normalize illumination.
            def inv_median(ccd_obj: CCDData) -> float:
                arr = np.ma.array(ccd_obj.data)
                med = np.ma.median(arr)
                if np.ma.is_masked(med) or float(med) == 0.0:
                    return 1.0
                return float(1.0 / med)

            combine_kwargs_flat = dict(
                method="median",
                scale=inv_median,
                sigma_clip=sigma,
                sigma_clip_func=np.ma.median,
                sigma_clip_dev_func=mad_std,
            )
            if mem_limit is not None:
                combine_kwargs_flat["mem_limit"] = float(mem_limit)

            master = combine(ccd_list, **combine_kwargs_flat)
            _add_history(
                master.meta,
                f"Master FLAT bin={binning} "
                f"gain={gain} exp={_float_str(exposure)}s temp={temperature}C",
            )
            _add_history(master.meta, f"Combined from {len(ccd_list)} files")

        if master is None:
            continue

        output.mkdir(parents=True, exist_ok=True)

        bin_str = _format_bin(binning)
        gain_str = f"g{gain}" if gain is not None else "gNA"
        exp_str = f"{_float_str(exposure)}s" if exposure is not None else "sNA"
        temp_str = f"{_float_str(temperature)}C" if temperature is not None else "CNA"

        out_name = output / f"master_{file_type_norm.lower()}_{bin_str}_{gain_str}_{exp_str}_{temp_str}.fits"

        try:
            master.data = master.data.astype(np.float32, copy=False)
            master.to_hdu(hdu_mask=None, hdu_uncertainty=None).writeto(
                out_name,
                overwrite=True,
            )
            logger.info(
                "Created master file %s from %d files", out_name, len(ccd_list)
            )
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
) -> List[Path]:
    """
    Calibrate frames from ``input_dir/<file_type>`` using given master frames.

    Calibration sequence
    --------------------
    For each group:
        1. Subtract master bias (if available and not no_bias).
        2. Subtract master dark (matched in exposure and temperature).
        3. Apply flat-field correction.

    Master selection
    ----------------
    Masters are discovered from ``master_files`` (file, directory,
    or iterable of paths), grouped with the same rules as raw frames.

    Parameters
    ----------
    input_dir : str or Path
        Root directory containing subdirectories by IMAGETYP.
    file_type : str
        Science-like image type to calibrate (e.g. 'LIGHT', 'OBJECT').
    output_dir : str or Path
        Output directory for calibrated FITS files.
    master_files : path-like or iterable
        Master calibration frames (or directory containing them).
    no_bias : bool
        If True, skip bias subtraction.
    dark_scale : bool
        If True, scale master darks by exposure time in subtract_dark.
    debug : bool
        If True, enable verbose logging.

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

    files = load_files(input_dir / file_type_norm, ("*.fit", "*.fits"))
    groups = grouping(files)

    master_paths = _normalize_master_paths(master_files)
    master_groups = grouping(master_paths)

    created: List[Path] = []

    for (imagetyp, binning, gain, exposure, temperature), file_list in groups.items():
        # Only treat frames whose header IMAGETYP matches requested file_type.
        if imagetyp != file_type_norm:
            continue
        if not file_list:
            continue

        # Find best master matches for this group.
        bias_match = None if no_bias else _find_single_match(
            master_groups, "BIAS", binning, gain, None, None
        )
        dark_match = _find_single_match(
            master_groups, "DARK", binning, gain, exposure, temperature
        )
        flat_match = _find_single_match(
            master_groups, "FLAT", binning, gain, None, None
        )

        master_bias = None
        master_dark = None
        master_flat = None

        if bias_match:
            bias_path, _ = bias_match
            with fits.open(bias_path) as hdul:
                master_bias = _to_ccd(hdul[0])

        if dark_match:
            dark_path, _ = dark_match
            with fits.open(dark_path) as hdul:
                master_dark = _to_ccd(hdul[0])

        if flat_match:
            flat_path, _ = flat_match
            with fits.open(flat_path) as hdul:
                master_flat = _to_ccd(hdul[0])

        for fp, _ in file_list:
            try:
                with fits.open(fp) as hdul:
                    ccd = _to_ccd(hdul[0])

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
                    # Set a floor to avoid division by zero / dead pixels
                    flat_data = np.asarray(master_flat.data, dtype=float)
                    finite = flat_data[np.isfinite(flat_data)]
                    if finite.size > 0:
                        # e.g. 1% din mediana valorilor finite
                        flat_floor = max(1e-6, 0.01 * float(np.median(finite)))
                    else:
                        flat_floor = 1e-6

                    ccd = flat_correct(ccd, master_flat, min_value=flat_floor)
                    _add_history(ccd.meta, f"Flat-field corrected (min_value={flat_floor:.3g})")
                else:
                    logger.warning("Missing FLAT for %s", fp)

                output_dir.mkdir(parents=True, exist_ok=True)

                out_path = output_dir / f"{Path(fp).stem}_cal.fits"
                ccd.data = ccd.data.astype(np.float32, copy=False)
                ccd.write(out_path, overwrite=True)
                created.append(out_path)

            except Exception:
                logger.exception("Calibration failed for %s", fp)

    return created
