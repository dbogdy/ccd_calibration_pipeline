#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Master DARK builder (v3, ccdproc-based).

Usage:
    python master_dark.py <path> [masters_folder] [no_bias]

With "no_bias", master bias lookup and subtraction are skipped entirely:
master darks are built from the raw dark frames (BIASSUB=False in the
header). Recommended for CMOS cameras with unstable bias levels - the
master dark then carries the bias inside it and calibrating a light with
a matching-exposure dark removes both.

Looks for a "Dark" folder (case-insensitive) inside <path>, groups the
DARK frames by binning, gain, exposure and CCD temperature, and writes one
master dark per group into <path>/<masters_folder> ("masters" default).

Processing per group (per the Astropy CCD Data Reduction Guide):
  1. Quality control  - level/noise outliers rejected before combining.
  2. Bias subtraction - ccdproc.subtract_bias with the matching master
     BIAS (same binning/gain and set-point temperature), guarded
     by a median-level screen: dark current is never negative, so when the
     group's predicted post-subtraction median is clearly negative (or most
     of the frames would be), the bias does not belong to this session and
     is NOT applied (recorded on screen and in the header). Individual
     frames dipping slightly below zero only produce a note - they are
     ordinary scatter around a level that is correct for the group.
  3. Combination      - ccdproc.combine: average + sigma clipping around
     the median (mad_std), memory-limited.
"""

import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

from calib_common import (collect_fits_files, env_float, fexp,
                          find_subfolder, tag)
from grouping import group_frames
from masters import find_master, index_masters
from ccd_ops import (MASTER_COMBINE_METHOD, MEM_BUDGET_MB, MIN_FRAMES_FOR_QC,
                     ReducedTempDir, apply_overscan_trim, combine_frames,
                     enough_frames, estimate_pair_noise,
                     filter_majority_shape, map_frames, print_group_stats,
                     print_qc_summary, quality_control, read_ccd,
                     scan_frames, warn_if_nonfinite)

import ccdproc as ccdp

# Dark current is never negative: when the group's predicted median after
# bias subtraction is below -this, the master bias does not match this dark
# session (CMOS bias levels drift/jump) and is not applied.
NEG_MEDIAN_TOL_ADU = env_float("CALIB_NEG_MEDIAN_TOL_ADU", 1.0)


def main():
    args = sys.argv[1:]
    skip_bias = False
    for flag in ("no_bias", "--no_bias", "--no-bias"):
        if flag in args:
            args.remove(flag)
            skip_bias = True

    if len(args) not in (1, 2):
        print(f"Usage: python {Path(sys.argv[0]).name} "
              "<path> [masters_folder] [no_bias]")
        sys.exit(1)

    root = Path(args[0])
    masters_name = args[1] if len(args) == 2 else "masters"
    if not root.is_dir():
        print(f"ERROR: path does not exist or is not a directory: {root}")
        sys.exit(1)

    dark_dir = find_subfolder(root, "Dark")
    if dark_dir is None:
        print(f"ERROR: no 'Dark' folder found in {root}")
        sys.exit(1)

    files = collect_fits_files(dark_dir)
    if not files:
        print(f"No .fit/.fits files found in {dark_dir}")
        sys.exit(1)

    print(f"Found {len(files)} FITS files in {dark_dir}")

    groups, skipped = group_frames(files, imagetyp="DARK",
                                   use_exposure=True, use_temp=True)
    if not groups:
        print("No DARK frames found - nothing to do.")
        sys.exit(1)

    n_frames = sum(len(g["files"]) for g in groups)
    print(f"{n_frames} DARK frames in {len(groups)} group(s)"
          f" ({skipped} file(s) skipped)")

    masters_dir = root / masters_name
    if skip_bias:
        bias_index = []
        print("Bias subtraction disabled (no_bias option) - "
              "master darks will be built from raw dark frames")
    else:
        bias_index = index_masters(masters_dir, "BIAS")
        if bias_index:
            print(f"Found {len(bias_index)} master BIAS file(s) in {masters_dir}")
        else:
            print(f"WARNING: no master BIAS files found in {masters_dir} - "
                  "master darks will be built without bias correction")
    masters_dir.mkdir(parents=True, exist_ok=True)

    sort_key = lambda g: (g["binx"] or 0, g["biny"] or 0, g["gain"] or 0,
                          g["exp"] if g["exp"] is not None else -1,
                          g["temp"] if g["temp"] is not None else 999)
    for group in sorted(groups, key=sort_key):
        binx, biny = group["binx"], group["biny"]
        gain, exp, temp = group["gain"], group["exp"], group["temp"]
        print(f"\nGroup bin={binx}x{biny} gain={gain} exp={fexp(exp)}s "
              f"temp={temp}C: {len(group['files'])} frames")

        # Cheap check first: a group too small to combine is not worth reading.
        if not enough_frames(len(group["files"]), group["files"]):
            continue

        frames = scan_frames(group["files"])
        if not frames:
            print("  No readable frames in this group, skipping.")
            continue
        frames, shape = filter_majority_shape(frames)
        print_group_stats(frames)

        # --- Quality control ---
        if len(frames) >= MIN_FRAMES_FOR_QC:
            frames, rejected, qc_info = quality_control(frames)
            print_qc_summary(len(frames), rejected, qc_info)
        else:
            print(f"  {'QC':<13}: skipped, only {len(frames)} frames "
                  f"(<{MIN_FRAMES_FOR_QC}) - no rejection applied")

        # ... and again once shape filtering and QC have taken their cut.
        if not enough_frames(len(frames), [fr["file"] for fr in frames],
                             stage="left after QC"):
            continue

        # --- Master bias selection + median-level damage screen ---
        bias_entry = None if skip_bias else find_master(
            bias_index, binx, biny, gain, temp, any_exposure=True)
        bias_ccd = None
        bias_note = None
        near_zero = None
        if skip_bias:
            bias_note = "bias subtraction disabled (no_bias option)"
            print(f"  {'Master BIAS':<13}: NOT APPLIED - {bias_note}")
        elif bias_entry is None:
            bias_note = "no matching master BIAS (binning/gain/temp)"
            print(f"  {'Master BIAS':<13}: NOT APPLIED - {bias_note}")
        else:
            bias_ccd = read_ccd(bias_entry["path"])
            reason = None
            if bias_ccd.data.shape != shape:
                reason = (f"shape mismatch: dark {shape} vs "
                          f"bias {bias_ccd.data.shape}")
            else:
                # Screen the group with the medians already collected:
                # median(dark) - median(bias) predicts the post-subtraction
                # level and costs no extra I/O. The verdict is taken on the
                # GROUP (its median, plus how many frames disagree), not on
                # the single worst frame: individual frames sitting one ADC
                # step low are ordinary scatter, not a mismatched bias.
                bias_median = float(np.median(bias_ccd.data))
                predicted = [fr["median"] - bias_median for fr in frames]
                group_predicted = float(np.median(predicted))
                bad = [p for p in predicted if p < -NEG_MEDIAN_TOL_ADU]
                if group_predicted < -NEG_MEDIAN_TOL_ADU \
                        or len(bad) > len(frames) / 2:
                    reason = (f"the group median would be "
                              f"{group_predicted:+.2f} ADU after subtraction "
                              f"and {len(bad)}/{len(frames)} frames would go "
                              f"below -{NEG_MEDIAN_TOL_ADU} (worst "
                              f"{min(predicted):.2f} ADU) - bias level does "
                              f"not match this dark session")
                elif bad:
                    near_zero = (f"{len(bad)}/{len(frames)} frames dip below "
                                 f"-{NEG_MEDIAN_TOL_ADU} ADU (worst "
                                 f"{min(predicted):.2f}) but the group median "
                                 f"stays at {group_predicted:+.2f} ADU - "
                                 f"bias applied")
            if reason:
                bias_note = (f"master BIAS {bias_entry['path'].name} "
                             f"would damage the darks: {reason}")
                print(f"  {'Master BIAS':<13}: NOT APPLIED - {bias_note}")
                bias_ccd = None
            else:
                print(f"  {'Master BIAS':<13}: {bias_entry['path'].name} "
                      f"(ccdproc.subtract_bias)")
                if near_zero:
                    print(f"  {'':<13}  NOTE: {near_zero}")

        use_files = [fr["file"] for fr in frames]
        first_header = fits.getheader(use_files[0], ext=0)
        needs_prep = (bias_ccd is not None
                      or bool(first_header.get("BIASSEC"))
                      or bool(first_header.get("TRIMSEC")))

        with ReducedTempDir(masters_dir, "dark") as tmp:
            if needs_prep:
                print(f"  {'Reducing':<13}: {len(use_files)} frames "
                      f"(ccdproc) before combining ...")

                def prep(f):
                    ccd, _ = apply_overscan_trim(read_ccd(f))
                    if bias_ccd is not None:
                        ccd = ccdp.subtract_bias(ccd, bias_ccd)
                    return tmp.write(ccd, f.name)

                use_files = map_frames(prep, use_files)

            pd_noise = estimate_pair_noise(use_files)
            if pd_noise is not None:
                print(f"  {'Pair noise':<13}: {pd_noise:.3f} ADU "
                      "(read + dark shot noise)")

            print(f"  {'Combining':<13}: {len(use_files)} frames, "
                  f"{MASTER_COMBINE_METHOD} + sigma clip, "
                  f"mem {MEM_BUDGET_MB} MB ...")
            master = combine_frames(use_files)

        master_median = float(np.median(master.data))

        master.header["IMAGETYP"] = "DARK"
        master.header["NCOMBINE"] = (len(use_files),
                                     "Number of frames combined")
        master.header["COMBINED"] = (True, "ccdproc combined master frame")
        master.header["DARKMED"] = (round(master_median, 4),
                                    "Median of master dark [ADU]")
        if master_median < -NEG_MEDIAN_TOL_ADU:
            warn = (f"master dark median is negative ({master_median:.2f} "
                    "ADU); dark current cannot be negative - the master "
                    "BIAS level likely does not match this dark session")
            print(f"  WARNING: {warn}")
            master.header.add_history(f"WARNING: {warn}")
        if pd_noise is not None:
            master.header["PDNOISE"] = (round(pd_noise, 4),
                                        "Pair-difference noise [ADU]")
        if bias_ccd is not None:
            master.header["BIASSUB"] = (True, "Master bias subtracted")
            master.header["BIASFILE"] = (bias_entry["path"].name,
                                         "Master bias used")
            master.header.add_history(
                f"Master bias subtracted (ccdproc.subtract_bias): "
                f"{bias_entry['path'].name}")
            if near_zero:
                master.header.add_history(f"NOTE: {near_zero}")
        else:
            master.header["BIASSUB"] = (False, "Master bias NOT subtracted")
            master.header.add_history(f"Master bias NOT applied: {bias_note}")
        master.header.add_history(
            f"Master DARK: ccdproc.combine {MASTER_COMBINE_METHOD} + sigma "
            f"clip of {len(use_files)} frames, bin={binx}x{biny} gain={gain} "
            f"exp={fexp(exp)}s temp={temp}C")
        warn_if_nonfinite(master)

        out_name = (
            f"master_dark_bin{tag(binx)}x{tag(biny)}"
            f"_gain{tag(gain)}_exp{fexp(exp)}s_temp{tag(temp)}C.fits"
        )
        out_path = masters_dir / out_name
        master.write(out_path, overwrite=True)
        print(f"  -> {out_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
