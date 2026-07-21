#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Master FLAT builder (v3, ccdproc-based).

Usage:
    python master_flat.py <path> [masters_folder] [dark_scale]

Looks for a "Flat" folder (case-insensitive) inside <path>, groups the
FLAT frames by binning, gain, exposure, CCD temperature and FILTER, and
writes one master flat per group into <path>/<masters_folder>.

Master selection per group (searched in the same masters folder):
  - master DARK: same binning/gain and set-point temperature, same
    exposure. If none matches the exposure and "dark_scale" is given,
    the closest-exposure dark is scaled by t_flat/t_dark - but ONLY when
    that dark is bias-subtracted (BIASSUB=True): bias does not scale with
    time. Otherwise the script stops with an error.
  - master BIAS: needed only when the chosen dark is bias-subtracted.
    When NO usable dark exists at all but a master BIAS matches, the bias
    alone is subtracted (an additive offset left in a flat dilutes its
    multiplicative structure).

Processing per group (per the Astropy CCD Data Reduction Guide):
  1. Each flat is reduced with ccdproc: overscan/trim (when present) ->
     subtract_bias -> subtract_dark (scale=True when allowed), written to
     a temporary folder.
  2. Quality control on the reduced frames: saturated / underexposed /
     non-positive medians are rejected (brightness drift is NOT - the
     per-frame scaling handles it).
  3. Combination: ccdproc.combine with scale=1/median (each frame counts
     equally), average + sigma clipping (mad_std), memory-limited. The
     result is renormalized to a median of exactly 1.0.
"""

import re
import sys
from pathlib import Path

import numpy as np
import astropy.units as u
import ccdproc as ccdp
from astropy.io import fits

from calib_common import (collect_fits_files, env_float, env_int, fexp,
                          find_subfolder, tag)
from grouping import group_frames
from masters import find_master, index_masters
from ccd_ops import (MASTER_COMBINE_METHOD, ReducedTempDir,
                     apply_overscan_trim, combine_frames,
                     filter_majority_shape, inv_median, print_frame_table,
                     read_ccd, scan_frames, warn_if_nonfinite)

# QC thresholds for flats (overridable via config.ini [master_flat])
MIN_FRAMES_FOR_QC = env_int("CALIB_FLAT_MIN_FRAMES_FOR_QC", 3)
# below this, only warn - never reject
SATURATION_ADU = env_float("CALIB_SATURATION_ADU", 65535.0)
# full scale of the camera ADC output
SAT_MEDIAN_FRAC = env_float("CALIB_SAT_MEDIAN_FRAC", 0.90)
# reject frames with median above this fraction
MIN_MEDIAN_FRAC = env_float("CALIB_MIN_MEDIAN_FRAC", 0.05)
# reject frames with median below this fraction
# (of saturation, measured after calibration)


def quality_control(frames):
    """
    Reject saturated, underexposed or over-corrected flats. Brightness
    drift between frames is normal for flats and is handled by the
    per-frame scaling, so median consistency is NOT checked.
    """
    accepted, rejected = [], []
    for f in frames:
        reason = None
        # non-positive first: it also falls below the underexposure cut,
        # but the real cause is over-correction, not exposure
        if f["median"] <= 0:
            reason = f"non-positive median ({f['median']:.2f}) after calibration"
        elif f["median"] > SAT_MEDIAN_FRAC * SATURATION_ADU:
            reason = (f"saturated: median {f['median']:.0f} > "
                      f"{SAT_MEDIAN_FRAC * 100:.0f}% of {SATURATION_ADU:.0f}")
        elif f["median"] < MIN_MEDIAN_FRAC * SATURATION_ADU:
            reason = (f"underexposed: median {f['median']:.0f} < "
                      f"{MIN_MEDIAN_FRAC * 100:.0f}% of {SATURATION_ADU:.0f}")
        if reason:
            rejected.append((f, reason))
        else:
            accepted.append(f)

    if not accepted:
        print("  WARNING: QC would reject every frame - keeping all")
        return frames, []
    return accepted, rejected


def main():
    args = sys.argv[1:]
    dark_scale = False
    for flag in ("dark_scale", "--dark_scale", "--dark-scale"):
        if flag in args:
            args.remove(flag)
            dark_scale = True

    if len(args) not in (1, 2):
        print(f"Usage: python {Path(sys.argv[0]).name} "
              "<path> [masters_folder] [dark_scale]")
        sys.exit(1)

    root = Path(args[0])
    masters_name = args[1] if len(args) == 2 else "masters"
    if not root.is_dir():
        print(f"ERROR: path does not exist or is not a directory: {root}")
        sys.exit(1)

    flat_dir = find_subfolder(root, "Flat")
    if flat_dir is None:
        print(f"ERROR: no 'Flat' folder found in {root}")
        sys.exit(1)

    files = collect_fits_files(flat_dir)
    if not files:
        print(f"No .fit/.fits files found in {flat_dir}")
        sys.exit(1)

    print(f"Found {len(files)} FITS files in {flat_dir}")

    groups, skipped = group_frames(files, imagetyp="FLAT", use_exposure=True,
                                   use_filter=True, use_temp=True)
    if not groups:
        print("No FLAT frames found - nothing to do.")
        sys.exit(1)

    n_frames = sum(len(g["files"]) for g in groups)
    print(f"{n_frames} FLAT frames in {len(groups)} group(s)"
          f" ({skipped} file(s) skipped)")

    masters_dir = root / masters_name
    dark_index = index_masters(masters_dir, "DARK")
    bias_index = index_masters(masters_dir, "BIAS")
    print(f"Found {len(dark_index)} master DARK and {len(bias_index)} "
          f"master BIAS file(s) in {masters_dir}")
    masters_dir.mkdir(parents=True, exist_ok=True)

    sort_key = lambda g: (g["binx"] or 0, g["biny"] or 0, g["gain"] or 0,
                          g["filter"] or "", g["exp"] if g["exp"] is not None else -1,
                          g["temp"] if g["temp"] is not None else 999)
    for group in sorted(groups, key=sort_key):
        binx, biny = group["binx"], group["biny"]
        gain, exp, temp, filt = group["gain"], group["exp"], group["temp"], group["filter"]
        filt_str = f" filter={filt}" if filt else ""
        print(f"\nGroup bin={binx}x{biny} gain={gain} exp={fexp(exp)}s "
              f"temp={temp}C{filt_str}: {len(group['files'])} frames")

        # --- Choose master dark (exact exposure, else scaled if allowed) ---
        dark_entry = find_master(dark_index, binx, biny, gain, temp, exp=exp)
        dark_scaled = False
        dark_note = None
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
                    dark_note = "no master DARK matches binning/gain/temp"
                    print(f"  NO DARK APPLIED: {dark_note}")
                else:
                    print(f"ERROR: dark scaling is NOT possible for group "
                          f"bin={binx}x{biny} gain={gain} exp={fexp(exp)}s "
                          f"temp={temp}C{filt_str}: no bias-subtracted master "
                          f"DARK is available (closest is "
                          f"{any_dark['path'].name}, {fexp(any_dark['exp'])}s, "
                          f"BIASSUB=False), and a dark can only be scaled to "
                          f"another exposure after the bias has been "
                          f"extracted from it. Rebuild the master darks with "
                          f"bias subtraction, or provide a master dark with "
                          f"exposure {fexp(exp)}s. Stopping.")
                    sys.exit(1)
            elif dark_entry["exp"] in (None, 0.0) or exp is None:
                dark_note = "cannot scale: missing exposure information"
                print(f"  NO DARK APPLIED: {dark_note}")
                dark_entry = None
            else:
                dark_scaled = True
                factor = exp / dark_entry["exp"]
                print(f"  Master DARK: {dark_entry['path'].name} "
                      f"scaled x{factor:.4f} "
                      f"({fexp(dark_entry['exp'])}s -> {fexp(exp)}s)")
        else:
            dark_note = (f"no master DARK with exposure {fexp(exp)}s "
                         "(run with dark_scale to allow scaling)")
            print(f"  NO DARK APPLIED: {dark_note}")

        # --- Bias: needed when the dark is bias-subtracted, or as the only
        # correction when no dark is available at all ---
        bias_entry = None
        if dark_entry is not None and dark_entry["biassub"]:
            bias_entry = find_master(bias_index, binx, biny, gain, temp,
                                     any_exposure=True)
            if bias_entry is None:
                print(f"ERROR: the master DARK {dark_entry['path'].name} is "
                      f"bias-subtracted (BIASSUB=True), so a master BIAS is "
                      f"required to calibrate the flats of group "
                      f"bin={binx}x{biny} gain={gain} exp={fexp(exp)}s "
                      f"temp={temp}C{filt_str} - but no master BIAS matching "
                      f"binning/gain/temperature was found in {masters_dir}. "
                      f"Stopping.")
                sys.exit(1)
            print(f"  Master BIAS: {bias_entry['path'].name}")
        elif dark_entry is not None:
            print("  Master BIAS not needed: the master DARK still contains "
                  "the bias (BIASSUB=False)")
        else:
            # No usable dark: subtract at least the bias when one matches.
            bias_entry = find_master(bias_index, binx, biny, gain, temp,
                                     any_exposure=True)
            if bias_entry is not None:
                print(f"  Master BIAS (dark unavailable): "
                      f"{bias_entry['path'].name} - subtracting bias only")
            else:
                print("  NO BIAS APPLIED: no matching master BIAS "
                      "(binning/gain/temp)")

        bias_ccd = read_ccd(bias_entry["path"]) if bias_entry else None
        dark_ccd = read_ccd(dark_entry["path"]) if dark_entry else None

        # --- Reduce each flat with ccdproc into a temp folder ---
        raw_frames = scan_frames(group["files"])
        if not raw_frames:
            print("  No readable frames in this group, skipping.")
            continue
        raw_frames, shape = filter_majority_shape(raw_frames)

        master_shapes = [c.data.shape for c in (bias_ccd, dark_ccd)
                         if c is not None]
        if any(s != shape for s in master_shapes):
            print(f"  WARNING: master shape does not match flats {shape} - "
                  "masters NOT applied")
            bias_ccd = dark_ccd = None
            bias_entry = dark_entry = None
            dark_scaled = False
            dark_note = (f"master shape mismatch vs flats {shape} - "
                         "bias/dark NOT applied")

        first_header = fits.getheader(raw_frames[0]["file"], ext=0)
        needs_prep = (bias_ccd is not None or dark_ccd is not None
                      or bool(first_header.get("BIASSEC"))
                      or bool(first_header.get("TRIMSEC")))

        with ReducedTempDir(masters_dir, "flat") as tmp:
            use_files = [fr["file"] for fr in raw_frames]
            if needs_prep:
                print("  Reducing each flat (ccdproc: overscan -> bias -> "
                      "dark) ...")
                prepped = []
                for f in use_files:
                    ccd = read_ccd(f)
                    ccd, _ = apply_overscan_trim(ccd)
                    if bias_ccd is not None:
                        ccd = ccdp.subtract_bias(ccd, bias_ccd)
                    if dark_ccd is not None:
                        ccd = ccdp.subtract_dark(
                            ccd, dark_ccd,
                            exposure_time="EXPTIME",
                            exposure_unit=u.second,
                            scale=dark_scaled)
                    prepped.append(tmp.write(ccd, f.name))
                use_files = prepped

            # --- QC on the frames as they will be combined ---
            frames = scan_frames(use_files)
            print_frame_table(frames)
            if len(frames) >= MIN_FRAMES_FOR_QC:
                frames, rejected = quality_control(frames)
                for fr, reason in rejected:
                    print(f"  REJECTED {fr['file'].name}: {reason}")
            elif len(frames) > 1:
                print(f"  QC: only {len(frames)} frames "
                      f"(<{MIN_FRAMES_FOR_QC}), no rejection applied")
            # scale=1/median cannot work on frames with median <= 0
            for fr in [fr for fr in frames if fr["median"] <= 0]:
                print(f"  REJECTED {fr['file'].name}: non-positive median "
                      f"({fr['median']:.2f}) - cannot scale")
            frames = [fr for fr in frames if fr["median"] > 0]
            if not frames:
                print("  Nothing left to combine, skipping group.")
                continue

            typical_level = float(np.median([fr["median"] for fr in frames]))

            print(f"  Combining {len(frames)} frames (ccdproc.combine, "
                  f"{MASTER_COMBINE_METHOD}, scale=1/median, "
                  f"typical level {typical_level:.0f} ADU) ...")
            master = combine_frames([fr["file"] for fr in frames],
                                    scale=inv_median)

        # Renormalize so the master's median is exactly 1.0.
        norm = float(np.median(master.data[np.isfinite(master.data)]))
        if norm > 0:
            master.data = master.data / norm
        master_min = float(np.nanmin(master.data))
        if master_min <= 0:
            n_bad = int(np.count_nonzero(master.data <= 0))
            print(f"  WARNING: master flat contains {n_bad} pixels <= 0 "
                  "(dead pixels or over-correction) - dividing by them "
                  "will misbehave; consider a bad-pixel map")

        master.header["IMAGETYP"] = "FLAT"
        master.header["NCOMBINE"] = (len(frames), "Number of frames combined")
        master.header["COMBINED"] = (True, "ccdproc combined master frame")
        master.header["FLATMED"] = (round(typical_level, 2),
                                    "Typical calibrated flat level [ADU]")
        master.header["DARKSUB"] = (dark_entry is not None,
                                    "Master dark subtracted")
        if dark_entry is not None:
            master.header["DARKFILE"] = (dark_entry["path"].name,
                                         "Master dark used")
            master.header["DARKSCAL"] = (dark_scaled,
                                         "Dark scaled by exposure ratio")
        master.header["BIASSUB"] = (bias_entry is not None,
                                    "Master bias subtracted from flats")
        if bias_entry is not None:
            master.header["BIASFILE"] = (bias_entry["path"].name,
                                         "Master bias used")
        master.header.add_history(
            f"Master FLAT: ccdproc.combine {MASTER_COMBINE_METHOD} + sigma "
            f"clip with scale=1/median of {len(frames)} frames, "
            f"bin={binx}x{biny} gain={gain} exp={fexp(exp)}s temp={temp}C"
            + (f" filter={filt}" if filt else ""))
        if dark_entry is not None:
            note = (f"Master dark subtracted (ccdproc.subtract_dark"
                    f"{', scaled' if dark_scaled else ''}): "
                    f"{dark_entry['path'].name}")
            master.header.add_history(note)
        if bias_entry is not None:
            master.header.add_history(
                f"Master bias subtracted (ccdproc.subtract_bias): "
                f"{bias_entry['path'].name}")
        if dark_note:
            master.header.add_history(f"NOTE: {dark_note}")
        master.header.add_history("Master flat normalized to median = 1.0")
        warn_if_nonfinite(master)

        filt_tag = ""
        if filt:
            safe = re.sub(r"[^A-Za-z0-9_+-]", "", filt)
            filt_tag = f"_filt{safe or 'NA'}"
        out_name = (
            f"master_flat_bin{tag(binx)}x{tag(biny)}"
            f"_gain{tag(gain)}_exp{fexp(exp)}s_temp{tag(temp)}C{filt_tag}.fits"
        )
        out_path = masters_dir / out_name
        master.write(out_path, overwrite=True)
        print(f"  -> {out_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
