#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Science-frame (LIGHT) calibration (v3, ccdproc-based).

Usage:
    python light_calibration.py <path> [masters_folder] [dark_scale]

Looks for a "Light" folder (case-insensitive) inside <path>, reads all
*.fit / *.fits science frames (IMAGETYP=LIGHT), groups them by binning,
gain, exposure, CCD temperature and FILTER, selects the required masters
from <path>/<masters_folder> ("masters" by default) and writes the
calibrated frames to <path>/calibrated, preserving the subfolder structure
(subfolders are typically object / panel names - for electroluminescence
work each panel usually has its own subfolder).

Calibration per frame, done entirely with ccdproc:

    ccdproc.subtract_overscan / trim_image     (when BIASSEC/TRIMSEC exist)
    ccdproc.subtract_bias                      (only when the dark is
                                                bias-subtracted)
    ccdproc.subtract_dark(scale=...)           (exact exposure, or scaled
                                                by the exposure ratio when
                                                "dark_scale" is allowed)
    ccdproc.flat_correct(min_value, norm_value)

Master selection rules (same as the master builders):
  - master DARK: same binning/gain and set-point temperature, same
    exposure. With "dark_scale", the closest-exposure BIAS-SUBTRACTED
    dark is scaled; scaling a dark that still contains bias is invalid and
    stops with an error. Missing dark = fatal error.
  - master BIAS: required (fatal if missing) only when the dark is
    bias-subtracted (BIASSUB=True).
  - master FLAT: same binning/gain/FILTER; temperature and exposure are
    not constraints (a normalized flat does not depend on them). The flat
    is renormalized to median 1.0 via flat_correct's norm_value; pixels
    below FLAT_MIN are clamped to FLAT_MIN by min_value (ccdproc's
    protection against division by dead/heavily vignetted pixels).
"""

import sys
from pathlib import Path

import numpy as np
import astropy.units as u
import ccdproc as ccdp

from calib_common import (as_float, collect_fits_files, env_float, env_path,
                          fexp, find_subfolder)
from grouping import group_frames
from masters import find_master, index_masters
from ccd_ops import (CR_OBJLIM, CR_READNOISE, CR_SIGCLIP, HOT_PIXEL_SIGMA,
                     apply_overscan_trim, compact, interpolate_bad_pixels,
                     make_hot_pixel_mask, read_ccd, reject_cosmic_rays)

# Flat handling (overridable via config.ini [calibration])
FLAT_MIN = env_float("CALIB_FLAT_MIN", 0.2)
# flat pixels below this are clamped (min_value)


def main():
    args = sys.argv[1:]
    dark_scale = False
    for flag in ("dark_scale", "--dark_scale", "--dark-scale"):
        if flag in args:
            args.remove(flag)
            dark_scale = True
    fix_hot_pixels = False
    for flag in ("fix_hot_pixels", "--fix_hot_pixels", "--fix-hot-pixels"):
        if flag in args:
            args.remove(flag)
            fix_hot_pixels = True
    fix_cosmic_rays = False
    for flag in ("fix_cosmic_rays", "--fix_cosmic_rays", "--fix-cosmic-rays"):
        if flag in args:
            args.remove(flag)
            fix_cosmic_rays = True

    if len(args) not in (1, 2):
        print(f"Usage: python {Path(sys.argv[0]).name} "
              "<path> [masters_folder] [dark_scale] "
              "[fix_hot_pixels] [fix_cosmic_rays]")
        sys.exit(1)

    root = Path(args[0])
    masters_name = args[1] if len(args) == 2 else "masters"
    if not root.is_dir():
        print(f"ERROR: path does not exist or is not a directory: {root}")
        sys.exit(1)

    light_dir = find_subfolder(root, "Light")
    if light_dir is None:
        print(f"ERROR: no 'Light' folder found in {root}")
        sys.exit(1)

    files = collect_fits_files(light_dir)
    if not files:
        print(f"No .fit/.fits files found in {light_dir}")
        sys.exit(1)

    print(f"Found {len(files)} FITS files in {light_dir}")

    groups, skipped = group_frames(files, imagetyp="LIGHT", use_exposure=True,
                                   use_filter=True, use_temp=True)
    if not groups:
        print("No LIGHT frames found - nothing to do.")
        sys.exit(1)

    n_frames = sum(len(g["files"]) for g in groups)
    print(f"{n_frames} LIGHT frames in {len(groups)} group(s)"
          f" ({skipped} file(s) skipped)")

    masters_dir = root / masters_name
    dark_index = index_masters(masters_dir, "DARK")
    bias_index = index_masters(masters_dir, "BIAS")
    flat_index = index_masters(masters_dir, "FLAT")
    print(f"Found {len(dark_index)} DARK, {len(bias_index)} BIAS, "
          f"{len(flat_index)} FLAT master(s) in {masters_dir}")

    # Output folder: CALIB_CALIBRATED_DIR (set by run_pipeline.py) or, when
    # run standalone, <path>/calibrated.
    out_root = env_path("CALIB_CALIBRATED_DIR", root / "calibrated")

    calibrated = 0
    sort_key = lambda g: (g["binx"] or 0, g["biny"] or 0, g["gain"] or 0,
                          g["filter"] or "", g["exp"] if g["exp"] is not None else -1,
                          g["temp"] if g["temp"] is not None else 999)
    for group in sorted(groups, key=sort_key):
        binx, biny = group["binx"], group["biny"]
        gain, exp, temp, filt = group["gain"], group["exp"], group["temp"], group["filter"]
        filt_str = f" filter={filt}" if filt else ""
        group_str = (f"bin={binx}x{biny} gain={gain} exp={fexp(exp)}s "
                     f"temp={temp}C{filt_str}")
        print(f"\nGroup {group_str}: {len(group['files'])} frames")

        # --- Master FLAT: binning/gain/filter; temp+exposure only break ties ---
        flat_entry = find_master(flat_index, binx, biny, gain, temp,
                                 any_exposure=True, filt=filt,
                                 match_filter=True, match_temp=False)
        if flat_entry is None:
            print(f"ERROR: no master FLAT matching binning/gain"
                  f"{'/filter ' + repr(filt) if filt else ''} was "
                  f"found in {masters_dir} for group {group_str}. A master "
                  f"flat is required to calibrate science frames. Stopping.")
            sys.exit(1)
        print(f"  Master FLAT: {flat_entry['path'].name}")

        # --- Master DARK: exact exposure, else scaled if allowed ---
        dark_entry = find_master(dark_index, binx, biny, gain, temp, exp=exp)
        dark_scaled = False
        if dark_entry is not None:
            print(f"  Master DARK: {dark_entry['path'].name}")
        elif dark_scale:
            # Only bias-subtracted darks can be scaled, so prefer those.
            scalable = [d for d in dark_index if d["biassub"]]
            dark_entry = find_master(scalable, binx, biny, gain, temp,
                                     exp=exp, any_exposure=True)
            if dark_entry is None:
                any_dark = find_master(dark_index, binx, biny, gain, temp,
                                       exp=exp, any_exposure=True)
                if any_dark is None:
                    print(f"ERROR: no master DARK matching binning/gain/"
                          f"temperature was found in {masters_dir} for group "
                          f"{group_str}. Stopping.")
                else:
                    print(f"ERROR: dark scaling is NOT possible for group "
                          f"{group_str}: no bias-subtracted master DARK is "
                          f"available (closest is {any_dark['path'].name}, "
                          f"{fexp(any_dark['exp'])}s, BIASSUB=False), and a "
                          f"dark can only be scaled to another exposure "
                          f"after the bias has been extracted from it. "
                          f"Rebuild the master darks with bias subtraction, "
                          f"or provide a master dark with exposure "
                          f"{fexp(exp)}s. Stopping.")
                sys.exit(1)
            if dark_entry["exp"] in (None, 0.0) or exp is None:
                print(f"ERROR: cannot scale master DARK for group "
                      f"{group_str}: missing exposure information. Stopping.")
                sys.exit(1)
            dark_scaled = True
            factor = exp / dark_entry["exp"]
            print(f"  Master DARK: {dark_entry['path'].name} "
                  f"scaled x{factor:.4f} "
                  f"({fexp(dark_entry['exp'])}s -> {fexp(exp)}s)")
        else:
            print(f"ERROR: no master DARK with exposure {fexp(exp)}s was "
                  f"found in {masters_dir} for group {group_str} "
                  f"(run with dark_scale to allow scaling a dark from "
                  f"another exposure). Stopping.")
            sys.exit(1)

        # --- Master BIAS: only when the dark is bias-subtracted ---
        bias_entry = None
        if dark_entry["biassub"]:
            bias_entry = find_master(bias_index, binx, biny, gain, temp,
                                     any_exposure=True)
            if bias_entry is None:
                print(f"ERROR: the master DARK {dark_entry['path'].name} is "
                      f"bias-subtracted (BIASSUB=True), so a master BIAS is "
                      f"required to calibrate group {group_str} - but no "
                      f"master BIAS matching binning/gain/temperature was "
                      f"found in {masters_dir}. Stopping.")
                sys.exit(1)
            print(f"  Master BIAS: {bias_entry['path'].name}")
        else:
            print("  Master BIAS not needed: the master DARK still contains "
                  "the bias (BIASSUB=False)")

        # --- Load masters once per group ---
        bias_ccd = read_ccd(bias_entry["path"]) if bias_entry else None
        dark_ccd = read_ccd(dark_entry["path"])
        flat_ccd = read_ccd(flat_entry["path"])

        flat_med = float(np.median(flat_ccd.data[np.isfinite(flat_ccd.data)]))
        if flat_med <= 0:
            print(f"ERROR: master flat {flat_entry['path'].name} has a "
                  f"non-positive median ({flat_med}) - unusable. Stopping.")
            sys.exit(1)
        if abs(flat_med - 1.0) > 0.01:
            print(f"  NOTE: master flat median is {flat_med:.4f}, "
                  "renormalizing to 1.0 during flat_correct")
        n_low = int(np.count_nonzero(
            ~np.isfinite(flat_ccd.data) | (flat_ccd.data < FLAT_MIN)))
        if n_low:
            print(f"  NOTE: {n_low} flat pixels < {FLAT_MIN} or non-finite "
                  f"will be clamped to {FLAT_MIN} (flat_correct min_value)")

        # --- Hot-pixel mask (optional): flagged once per group from the
        # master dark, then interpolated in every calibrated frame ---
        hot_mask = None
        if fix_hot_pixels:
            hot_mask = make_hot_pixel_mask(dark_ccd.data, sigma=HOT_PIXEL_SIGMA)
            print(f"  Hot-pixel mask: {int(hot_mask.sum())} pixel(s) "
                  f"(sigma={HOT_PIXEL_SIGMA}) from {dark_entry['path'].name}")

        # --- Calibrate each frame with ccdproc ---
        for f in group["files"]:
            try:
                ccd = read_ccd(f)
            except Exception as exc:
                print(f"  WARNING: cannot read {f.name}: {exc} - skipped")
                continue

            if ccd.data.shape != dark_ccd.data.shape \
                    or ccd.data.shape != flat_ccd.data.shape:
                print(f"  WARNING: {f.name} shape {ccd.data.shape} does not "
                      f"match masters - skipped")
                continue

            ccd, overscan_used = apply_overscan_trim(ccd)
            med_before = float(np.median(ccd.data))

            if bias_ccd is not None:
                ccd = ccdp.subtract_bias(ccd, bias_ccd)
            ccd = ccdp.subtract_dark(ccd, dark_ccd,
                                     exposure_time="EXPTIME",
                                     exposure_unit=u.second,
                                     scale=dark_scaled)
            ccd = ccdp.flat_correct(ccd, flat_ccd,
                                    min_value=FLAT_MIN,
                                    norm_value=flat_med)

            # --- Optional cleaning: hot-pixel repair, then cosmic rays ---
            if hot_mask is not None and hot_mask.any():
                ccd.data = interpolate_bad_pixels(ccd.data, hot_mask)
                ccd.header.add_history(
                    f"Hot pixel correction: {int(hot_mask.sum())} pixel(s) "
                    f"interpolated (from master dark)")
            n_cr = 0
            if fix_cosmic_rays:
                gain_val = max(0.1, as_float(ccd.header.get("EGAIN"))
                               or as_float(ccd.header.get("GAIN")) or 1.0)
                rn_val = max(0.1, as_float(ccd.header.get("RDNOISE"))
                             or CR_READNOISE)
                try:
                    ccd, n_cr = reject_cosmic_rays(ccd, gain_val, rn_val,
                                                   CR_SIGCLIP, CR_OBJLIM)
                    if n_cr:
                        ccd.header.add_history(
                            f"Cosmic ray rejection: {n_cr} pixel(s) cleaned "
                            f"(L.A.Cosmic)")
                except Exception as exc:
                    print(f"  WARNING: cosmic ray rejection failed for "
                          f"{f.name}: {exc}")

            med_after = float(np.median(ccd.data))

            ccd.header["BIASSUB"] = (bias_entry is not None,
                                     "Master bias subtracted")
            if bias_entry is not None:
                ccd.header["BIASFILE"] = (bias_entry["path"].name,
                                          "Master bias used")
            ccd.header["DARKFILE"] = (dark_entry["path"].name,
                                      "Master dark used")
            ccd.header["DARKSCAL"] = (dark_scaled,
                                      "Dark scaled by exposure ratio")
            ccd.header["FLATFILE"] = (flat_entry["path"].name,
                                      "Master flat used")
            ccd.header.add_history(
                f"Calibrated with ccdproc: "
                f"{'subtract_bias -> ' if bias_entry else ''}"
                f"subtract_dark{'(scaled)' if dark_scaled else ''} -> "
                f"flat_correct(min_value={FLAT_MIN}, "
                f"norm_value={flat_med:.4f})")
            if bias_entry is not None:
                ccd.header.add_history(
                    f"Master bias: {bias_entry['path'].name}")
            ccd.header.add_history(f"Master dark: {dark_entry['path'].name}")
            ccd.header.add_history(f"Master flat: {flat_entry['path'].name}")
            if overscan_used:
                ccd.header.add_history(
                    "Overscan (BIASSEC) subtracted via ccdproc")

            rel = f.relative_to(light_dir)
            out_path = out_root / rel.parent / f"{f.stem}_cal.fits"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            compact(ccd).write(out_path, overwrite=True)
            calibrated += 1
            cr_str = f" [CR:{n_cr}]" if fix_cosmic_rays and n_cr else ""
            print(f"  {f.name}: median {med_before:.1f} -> {med_after:.1f}"
                  f"{cr_str}  -> {out_path.relative_to(root)}")

    print(f"\nDone. {calibrated} frame(s) calibrated into {out_root}")


if __name__ == "__main__":
    main()
