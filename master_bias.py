#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Master BIAS builder (v3, ccdproc-based).

Usage:
    python master_bias.py <path> [output_folder]

Looks for a "Bias" folder (case-insensitive) inside <path>, loads all
*.fit / *.fits files whose IMAGETYP is BIAS, groups them by binning, gain
and set-point temperature (SET-TEMP; frames whose CCD-TEMP deviates more
than TEMP_MATCH_TOL from it are dropped), and writes one master bias per
group into <path>/<output_folder> ("masters" by default).

Processing per group (per the Astropy CCD Data Reduction Guide):
  1. Overscan/trim   - ccdproc.subtract_overscan / trim_image when the
     headers carry BIASSEC / TRIMSEC (skipped otherwise; in that case the
     original files are combined directly, without temporary copies).
  2. Quality control - frames whose median level or noise deviates strongly
     from the rest of the group are rejected before combining.
  3. Read noise      - pair-difference estimate robust_std(a-b)/sqrt(2),
     stored as RDNOISE [ADU] (and RDNOISEE [e-] when EGAIN is present).
  4. Combination     - ccdproc.combine: average with sigma clipping around
     the median (mad_std deviation) and a memory limit, i.e. the guide's
     recommended master-bias recipe.
"""

import sys
from pathlib import Path

from astropy.io import fits

from calib_common import collect_fits_files, find_subfolder, tag
from grouping import group_frames
from ccd_ops import (MASTER_COMBINE_METHOD, MIN_FRAMES_FOR_QC,
                     ReducedTempDir, apply_overscan_trim, combine_frames,
                     estimate_pair_noise, filter_majority_shape,
                     print_frame_table, quality_control, read_ccd,
                     scan_frames, warn_if_nonfinite)


def main():
    if len(sys.argv) not in (2, 3):
        print(f"Usage: python {Path(sys.argv[0]).name} <path> [output_folder]")
        sys.exit(1)

    root = Path(sys.argv[1])
    output_name = sys.argv[2] if len(sys.argv) == 3 else "masters"
    if not root.is_dir():
        print(f"ERROR: path does not exist or is not a directory: {root}")
        sys.exit(1)

    bias_dir = find_subfolder(root, "Bias")
    if bias_dir is None:
        print(f"ERROR: no 'Bias' folder found in {root}")
        sys.exit(1)

    files = collect_fits_files(bias_dir)
    if not files:
        print(f"No .fit/.fits files found in {bias_dir}")
        sys.exit(1)

    print(f"Found {len(files)} FITS files in {bias_dir}")

    groups, skipped = group_frames(files, imagetyp="BIAS", use_temp=True)
    if not groups:
        print("No BIAS frames found - nothing to do.")
        sys.exit(1)

    n_frames = sum(len(g["files"]) for g in groups)
    print(f"{n_frames} BIAS frames in {len(groups)} group(s)"
          f" ({skipped} file(s) skipped)")

    masters_dir = root / output_name
    masters_dir.mkdir(parents=True, exist_ok=True)

    sort_key = lambda g: (g["binx"] or 0, g["biny"] or 0,
                          g["gain"] or 0, g["temp"] if g["temp"] is not None else 999)
    for group in sorted(groups, key=sort_key):
        binx, biny = group["binx"], group["biny"]
        gain, temp = group["gain"], group["temp"]
        print(f"\nGroup bin={binx}x{biny} gain={gain} temp={temp}C: "
              f"{len(group['files'])} frames")

        frames = scan_frames(group["files"])
        if not frames:
            print("  No readable frames in this group, skipping.")
            continue
        frames, shape = filter_majority_shape(frames)
        print_frame_table(frames)

        if len(frames) >= MIN_FRAMES_FOR_QC:
            frames, rejected = quality_control(frames)
            for fr, reason in rejected:
                print(f"  REJECTED {fr['file'].name}: {reason}")
        elif len(frames) > 1:
            print(f"  QC: only {len(frames)} frames (<{MIN_FRAMES_FOR_QC}), "
                  "no rejection applied")

        use_files = [fr["file"] for fr in frames]
        first_header = fits.getheader(use_files[0], ext=0)

        # Overscan/trim only when the camera actually provides the sections.
        needs_prep = bool(first_header.get("BIASSEC")
                          or first_header.get("TRIMSEC"))

        with ReducedTempDir(masters_dir, "bias") as tmp:
            if needs_prep:
                print("  Applying overscan/trim (ccdproc) to each frame ...")
                prepped = []
                for f in use_files:
                    ccd = read_ccd(f)
                    ccd, _ = apply_overscan_trim(ccd)
                    prepped.append(tmp.write(ccd, f.name))
                use_files = prepped

            # Read noise from the frames actually being combined.
            rn_adu = estimate_pair_noise(use_files)
            rn_e = None
            if rn_adu is not None:
                try:
                    rn_e = rn_adu * float(first_header.get("EGAIN"))
                except (TypeError, ValueError):
                    pass
                msg = f"  Read noise: {rn_adu:.3f} ADU"
                if rn_e is not None:
                    msg += f" = {rn_e:.3f} e- (via EGAIN)"
                print(msg)

            print(f"  Combining {len(use_files)} frames (ccdproc.combine, "
                  f"{MASTER_COMBINE_METHOD} + sigma clip) ...")
            master = combine_frames(use_files)

        master.header["IMAGETYP"] = "BIAS"
        master.header["NCOMBINE"] = (len(use_files),
                                     "Number of frames combined")
        master.header["COMBINED"] = (True, "ccdproc combined master frame")
        if rn_adu is not None:
            master.header["RDNOISE"] = (round(rn_adu, 4),
                                        "Read noise estimate [ADU]")
        if rn_e is not None:
            master.header["RDNOISEE"] = (round(rn_e, 4),
                                         "Read noise estimate [e-]")
        master.header.add_history(
            f"Master BIAS: ccdproc.combine {MASTER_COMBINE_METHOD} + sigma "
            f"clip of {len(use_files)} frames, bin={binx}x{biny} gain={gain} "
            f"temp={temp}C")
        if needs_prep:
            master.header.add_history(
                "Per-frame overscan/trim applied (ccdproc)")
        warn_if_nonfinite(master)

        out_name = (
            f"master_bias_bin{tag(binx)}x{tag(biny)}"
            f"_gain{tag(gain)}_temp{tag(temp)}C.fits"
        )
        out_path = masters_dir / out_name
        master.write(out_path, overwrite=True)
        print(f"  -> {out_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
