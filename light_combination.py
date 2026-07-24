#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combination of calibrated science frames (v3, ccdproc-based).

Usage:
    python light_combination.py <calibrated_path> [method]

Scans <calibrated_path> (the "calibrated" folder produced by
light_calibration.py) recursively, groups the frames and combines every
group that has more than one image. Results go to a "combined" folder
created NEXT TO the calibrated folder, preserving the subfolder structure
(subfolders are object/panel names and are never mixed).

Grouping key: subfolder, binning, gain, exposure, set-point temperature
(SET-TEMP, exact), FILTER and EL injection current.

FILTER / current handling (electroluminescence sessions):
  - Both are read from the FITS header first (FILTER and ELCURR keywords).
  - When absent there, they are parsed from the filename (tokens like
    "_5A_" for the current and the token that follows it for the filter,
    e.g. Light_71102_10.0s_..._-14.8C_5A_SDSS-y_0012_cal.fits).
  - Values parsed from the filename are written back into the calibrated
    file's header (only when the keyword is missing) and into the combined
    image's header. On disagreement the header wins (with a warning).

Combination methods (second argument, default "mean"), all executed with
ccdproc.combine (memory-limited, chunked I/O):
    mean    - method='average'
    median  - method='median' ("medium" accepted as alias)
    sum     - method='sum' (EXPTIME of the result = total exposure)
    sclip   - method='average' with 3-sigma clipping around the median
              (mad_std deviation) - best outlier rejection for >= ~5 frames
"""

import re
import sys
from pathlib import Path

from astropy.io import fits

from calib_common import collect_fits_files, env_path, fexp, tag
from grouping import group_frames
from ccd_ops import combine_frames

METHODS = ("mean", "median", "sum", "sclip")
ALIASES = {"medium": "median", "avg": "mean", "average": "mean"}


def parse_filename_tokens(name: str):
    """
    Extract (current, filter) from a calibrated filename.

    The current is a token like "5A" or "11A"; the filter is the token that
    immediately follows it (when it is not purely numeric), e.g.:
        Light_71102_10.0s_..._-14.8C_5A_SDSS-y_0012_cal.fits
        -> current "5A", filter "SDSS-y"
    Returns (None, None) parts when not found.
    """
    stem = Path(name).stem
    tokens = stem.split("_")
    current = None
    filt = None
    for i, tok in enumerate(tokens):
        if re.fullmatch(r"\d+(?:\.\d+)?A", tok):
            current = tok
            if i + 1 < len(tokens):
                nxt = tokens[i + 1]
                if nxt.lower() not in ("cal",) and not re.fullmatch(r"\d+", nxt):
                    filt = nxt
            break
    return current, filt


def ensure_header_keywords(path: Path, header):
    """
    Fill FILTER / ELCURR in a calibrated file's header from its filename
    when the keywords are missing. Header values always take priority; a
    mismatch only prints a warning. Returns (filter, current) - the
    effective values after reconciliation.
    """
    name_curr, name_filt = parse_filename_tokens(path.name)
    hdr_filt = str(header.get("FILTER", "")).strip() or None
    hdr_curr = str(header.get("ELCURR", "")).strip() or None

    updates = {}
    if hdr_filt is None and name_filt:
        updates["FILTER"] = (name_filt, "Filter (recovered from filename)")
        hdr_filt = name_filt
    elif hdr_filt and name_filt and hdr_filt != name_filt:
        print(f"  WARNING: {path.name}: header FILTER='{hdr_filt}' differs "
              f"from filename '{name_filt}' - using the header value")

    if hdr_curr is None and name_curr:
        updates["ELCURR"] = (name_curr, "EL injection current (from filename)")
        hdr_curr = name_curr
    elif hdr_curr and name_curr and hdr_curr != name_curr:
        print(f"  WARNING: {path.name}: header ELCURR='{hdr_curr}' differs "
              f"from filename '{name_curr}' - using the header value")

    if updates:
        try:
            with fits.open(path, mode="update") as hdul:
                for key, val in updates.items():
                    hdul[0].header[key] = val
        except Exception as exc:
            print(f"  WARNING: cannot update header of {path.name}: {exc}")

    return hdr_filt, hdr_curr


def main():
    args = sys.argv[1:]
    if len(args) not in (1, 2):
        print(f"Usage: python {Path(sys.argv[0]).name} "
              f"<calibrated_path> [method: {'/'.join(METHODS)}]")
        sys.exit(1)

    cal_dir = Path(args[0])
    method = args[1].lower() if len(args) == 2 else "mean"
    method = ALIASES.get(method, method)
    if method not in METHODS:
        print(f"ERROR: unknown combination method '{args[1]}'. "
              f"Available: {', '.join(METHODS)} (aliases: "
              f"{', '.join(f'{a}->{b}' for a, b in ALIASES.items())})")
        sys.exit(1)
    if not cal_dir.is_dir():
        print(f"ERROR: path does not exist or is not a directory: {cal_dir}")
        sys.exit(1)

    files = collect_fits_files(cal_dir)
    if not files:
        print(f"No .fit/.fits files found in {cal_dir}")
        sys.exit(1)

    print(f"Found {len(files)} FITS files in {cal_dir}")
    print(f"Combination method: {method}")

    # Extra exact-match keys for the generic grouping: subfolder, FILTER and
    # EL current (the last two reconciled between header and filename).
    def extra_keys(path, header):
        filt, curr = ensure_header_keywords(path, header)
        return {"subdir": path.parent.relative_to(cal_dir),
                "filter": filt, "current": curr}

    groups, _ = group_frames(files, use_exposure=True, use_temp=True,
                             extra_key_func=extra_keys)
    multi = [g for g in groups if len(g["files"]) > 1]
    single = [g for g in groups if len(g["files"]) == 1]
    print(f"{len(files)} frames in {len(groups)} group(s): "
          f"{len(multi)} to combine, {len(single)} single-frame (skipped)")
    for g in single:
        print(f"  single frame, not combined: {g['files'][0].name}")

    # Output folder: CALIB_COMBINED_DIR (set by run_pipeline.py) or, when run
    # standalone, a "combined" folder next to the calibrated one.
    out_root = env_path("CALIB_COMBINED_DIR", cal_dir.parent / "combined")

    combined = 0
    sort_key = lambda g: (str(g["subdir"]), g["filter"] or "",
                          g["current"] or "", g["exp"] if g["exp"] is not None else -1,
                          g["temp"] if g["temp"] is not None else 999)
    for g in sorted(multi, key=sort_key):
        parts = [f"bin={g['binx']}x{g['biny']}", f"gain={g['gain']}",
                 f"exp={fexp(g['exp'])}s", f"temp={g['temp']}C"]
        if g["filter"]:
            parts.append(f"filter={g['filter']}")
        if g["current"]:
            parts.append(f"current={g['current']}")
        sub = f"[{g['subdir']}] " if str(g["subdir"]) != "." else ""
        print(f"\nGroup {sub}{' '.join(parts)}: {len(g['files'])} frames")

        # All frames in a group must share one shape; drop the odd ones out.
        shapes = {}
        for f in g["files"]:
            try:
                shp = tuple(fits.getheader(f, ext=0)[k]
                            for k in ("NAXIS2", "NAXIS1"))
            except Exception:
                shp = None
            shapes.setdefault(shp, []).append(f)
        shape = max(shapes, key=lambda s: len(shapes[s]))
        use_files = shapes[shape]
        for s, fl in shapes.items():
            if s != shape:
                for f in fl:
                    print(f"  WARNING: {f.name} has shape {s}, "
                          f"group uses {shape} - excluded")
        if len(use_files) < 2:
            print("  Fewer than 2 usable frames left, skipping.")
            continue

        if method == "sclip":
            result = combine_frames(use_files, method="average",
                                    sigma_clip=True,
                                    sigma_low=3.0, sigma_high=3.0)
        else:
            ccd_method = {"mean": "average", "median": "median",
                          "sum": "sum"}[method]
            result = combine_frames(use_files, method=ccd_method,
                                    sigma_clip=False)

        header = result.header
        header["NCOMBINE"] = (len(use_files), "Number of frames combined")
        header["COMBMETH"] = (method, "Combination method")
        total_exp = (g["exp"] or 0) * len(use_files)
        header["EXPTOTAL"] = (round(total_exp, 6),
                              "Total exposure of combined frames [s]")
        if method == "sum" and g["exp"] is not None:
            header["EXPTIME"] = (total_exp,
                                 "Total exposure (sum combination) [s]")
        if g["filter"] and not str(header.get("FILTER", "")).strip():
            header["FILTER"] = (g["filter"], "Filter (recovered from filename)")
        if g["current"] and not str(header.get("ELCURR", "")).strip():
            header["ELCURR"] = (g["current"],
                                "EL injection current (from filename)")
        header.add_history(
            f"Combined {len(use_files)} calibrated frames with "
            f"ccdproc.combine, method={method}")
        header.add_history(f"Group: {' '.join(parts)}")

        filt_tag = ""
        if g["filter"]:
            safe = re.sub(r"[^A-Za-z0-9_+-]", "", g["filter"])
            filt_tag = f"_{safe}"
        curr_tag = f"_{g['current']}" if g["current"] else ""
        out_name = (
            f"combined_{method}_bin{tag(g['binx'])}x{tag(g['biny'])}"
            f"_gain{tag(g['gain'])}_exp{fexp(g['exp'])}s"
            f"_temp{tag(g['temp'])}C{curr_tag}{filt_tag}.fits"
        )
        out_path = out_root / g["subdir"] / out_name
        out_path.parent.mkdir(parents=True, exist_ok=True)
        result.write(out_path, overwrite=True)
        combined += 1
        print(f"  -> {out_path}")

    print(f"\nDone. {combined} combined image(s) written to {out_root}")


if __name__ == "__main__":
    main()
