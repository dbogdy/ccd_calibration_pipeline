#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Master-frame lookup shared by the pipeline scripts.

index_masters scans a masters folder and reads the matching keys from each
master's header; find_master picks the best candidate for a group of frames
(binning/gain exact, same set-point temperature, optional exposure and
filter constraints).

Temperature: like the grouping, masters are keyed on SET-TEMP (the
commanded temperature, kept in the master's header from its first frame),
falling back to the rounded CCD-TEMP for older masters. Both the group and
the master carry a discrete set point, so the match is exact - no
tolerance is needed.
"""

from pathlib import Path

from astropy.io import fits

from calib_common import as_int, as_float, collect_fits_files, exposures_match


def index_masters(masters_dir: Path, imagetyp: str):
    """
    Scan the masters folder for master files of the given IMAGETYP
    (e.g. "BIAS", "DARK", "FLAT").

    Returns a list of dicts: path, binx, biny, gain, exp, temp, filter,
    biassub. temp comes from SET-TEMP (fallback: rounded CCD-TEMP).
    """
    if not masters_dir.is_dir():
        return []
    index = []
    for f in collect_fits_files(masters_dir):
        try:
            header = fits.getheader(f, ext=0)
        except Exception:
            continue
        if str(header.get("IMAGETYP", "")).strip().upper() != imagetyp:
            continue
        set_t = as_int(header.get("SET-TEMP"))
        index.append({
            "path": f,
            "binx": as_int(header.get("XBINNING")),
            "biny": as_int(header.get("YBINNING")),
            "gain": as_int(header.get("GAIN")),
            "exp": as_float(header.get("EXPTIME")),
            "temp": set_t if set_t is not None
            else as_int(header.get("CCD-TEMP")),
            "filter": str(header.get("FILTER", "")).strip() or None,
            "biassub": bool(header.get("BIASSUB", False)),
        })
    return index


def find_master(index, binx, biny, gain, temp, exp=None, any_exposure=False,
                filt=None, match_filter=False, match_temp=True):
    """
    Pick the master matching binning and gain exactly and - with
    match_temp=True - the same set-point temperature (None matches None).
    With exp given, only exact exposure matches qualify; with
    any_exposure=True the exposure is unconstrained and the closest one
    wins. With match_filter=True the FILTER must match exactly (None
    matches None). With match_temp=False the temperature is not a
    constraint, only a tie-breaker (used for master flats: already
    calibrated and normalized, hence temperature-independent).
    Ties: closest temperature, then newest file. Returns the index entry or
    None.
    """
    candidates = []
    for b in index:
        if (b["binx"], b["biny"], b["gain"]) != (binx, biny, gain):
            continue
        if match_filter and b.get("filter") != filt:
            continue
        if match_temp:
            if b["temp"] != temp:
                continue
            dt = 0
        elif b["temp"] is not None and temp is not None:
            dt = abs(b["temp"] - temp)
        else:
            dt = 999
        if not any_exposure:
            if not exposures_match(b["exp"], exp):
                continue
            de = 0.0
        else:
            if b["exp"] is None or exp is None:
                de = 1e9
            else:
                de = abs(b["exp"] - exp)
        candidates.append((de, dt, b))
    if not candidates:
        return None
    best_de, best_dt, _ = min(candidates, key=lambda c: (c[0], c[1]))
    best = [b for de, dt, b in candidates if de == best_de and dt == best_dt]
    return max(best, key=lambda b: b["path"].stat().st_mtime)
