#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Frame grouping shared by every pipeline script.

All scripts group frames by (binning, gain) and, by default, by the
commanded sensor temperature; darks/flats/lights additionally match the
exposure, flats/lights the FILTER, and the combination script adds its own
exact-match keys (subfolder, EL current) through extra_key_func.

Temperature handling: groups are keyed on SET-TEMP (the temperature the
cooler was commanded to), not on the instantaneous CCD-TEMP - the set
point is a discrete value, so all frames of one set match exactly and no
tolerance chaining is needed. CCD-TEMP is used as a consistency check
instead: a frame whose |CCD-TEMP - SET-TEMP| > TEMP_MATCH_TOL was taken
before the cooler stabilized and is dropped entirely (with a message).
Frames without SET-TEMP fall back to the rounded CCD-TEMP as their group
temperature (no check possible).
"""

from astropy.io import fits

from calib_common import as_int, as_float, env_float, exposures_match

# deg C: max |CCD-TEMP - SET-TEMP| for a usable frame; overridable from
# config.ini ([general] temp_match_tol) via run_pipeline.py.
TEMP_MATCH_TOL = env_float("CALIB_TEMP_MATCH_TOL", 1.0)


def group_frames(files, imagetyp=None, use_exposure=False, use_filter=False,
                 use_temp=False, extra_key_func=None,
                 temp_match_tol=TEMP_MATCH_TOL):
    """
    Group frames by (binning, gain[, temperature][, exposure][, filter]
    [, extra keys]); all keys are matched exactly.

    imagetyp       - when given, frames whose IMAGETYP differs are skipped
                     (with a message) and counted in `skipped`.
    use_exposure   - match EXPTIME too (exact, via exposures_match).
    use_filter     - match the FILTER header keyword too (missing -> None).
    use_temp       - group by SET-TEMP (fallback: rounded CCD-TEMP) and drop
                     frames with |CCD-TEMP - SET-TEMP| > temp_match_tol
                     (cooler not stabilized); they are counted in `skipped`.
    temp_match_tol - threshold [deg C] for that check (default 1.0, or the
                     config.ini value when run through run_pipeline.py).
    extra_key_func - optional callable (path, header) -> dict of additional
                     exact-match keys; they are stored in the group dict.

    Returns (groups, skipped) - groups is a list of dicts with keys
    binx, biny, gain, temp, files (+ exp, filter and the extra keys when
    requested).
    """
    groups = []
    skipped = 0
    for f in files:
        try:
            header = fits.getheader(f, ext=0)
        except Exception as exc:
            print(f"  WARNING: cannot read header of {f.name}: {exc}")
            skipped += 1
            continue

        if imagetyp is not None:
            actual = str(header.get("IMAGETYP", "")).strip().upper()
            if actual != imagetyp:
                print(f"  Skipping {f.name}: IMAGETYP is "
                      f"'{actual or 'missing'}', not {imagetyp}")
                skipped += 1
                continue

        temp = None
        if use_temp:
            set_t = as_float(header.get("SET-TEMP"))
            ccd_t = as_float(header.get("CCD-TEMP"))
            if set_t is not None and ccd_t is not None \
                    and abs(ccd_t - set_t) > temp_match_tol:
                print(f"  Skipping {f.name}: CCD-TEMP {ccd_t:.1f}C deviates "
                      f"{abs(ccd_t - set_t):.1f}C from SET-TEMP {set_t:.1f}C "
                      f"(> {temp_match_tol}C - cooler not stabilized)")
                skipped += 1
                continue
            # SET-TEMP is the group key; frames without it fall back to the
            # rounded instantaneous temperature.
            temp = as_int(set_t) if set_t is not None else as_int(ccd_t)

        binx = as_int(header.get("XBINNING"))
        biny = as_int(header.get("YBINNING"))
        gain = as_int(header.get("GAIN"))
        exp = as_float(header.get("EXPTIME")) if use_exposure else None
        filt = (str(header.get("FILTER", "")).strip() or None) \
            if use_filter else None
        extra = extra_key_func(f, header) if extra_key_func else {}

        target = None
        for g in groups:
            if (g["binx"], g["biny"], g["gain"]) != (binx, biny, gain):
                continue
            if use_temp and g["temp"] != temp:
                continue
            if use_filter and g["filter"] != filt:
                continue
            if any(g[k] != v for k, v in extra.items()):
                continue
            if use_exposure and not exposures_match(g["exp"], exp):
                continue
            target = g
            break
        if target is None:
            target = {"binx": binx, "biny": biny, "gain": gain,
                      "temp": temp, "files": []}
            if use_exposure:
                target["exp"] = exp
            if use_filter:
                target["filter"] = filt
            target.update(extra)
            groups.append(target)
        target["files"].append(f)

    return groups, skipped
