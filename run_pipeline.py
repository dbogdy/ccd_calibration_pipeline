#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automatic runner for the calibration pipeline, driven by config.ini.

Usage:
    python run_pipeline.py [config_file]

Reads the configuration ("config.ini" next to this script by default),
then runs the enabled steps in order by invoking the existing scripts:

    1. master_bias.py        <images_path> <masters_dir>
    2. master_dark.py        <images_path> <masters_dir> [no_bias] [level_align]
    3. master_flat.py        <images_path> <masters_dir> [dark_scale]
    4. light_calibration.py  <images_path> <masters_dir> [dark_scale]
    5. light_combination.py  <calibrated_dir> [method]

Output layout: everything is written under one output folder. Set
[paths] output_folder (relative to images_path, or absolute) and the
masters, calibrated and combined frames land in subfolders of it
(<output_folder>/masters, /calibrated, /combined by default). The
calibrated/combined locations are passed to the scripts through the
CALIB_CALIBRATED_DIR / CALIB_COMBINED_DIR environment variables.

Each step inherits this process's console, so all output is shown live.
The pipeline stops at the first step that exits with an error, so a later
step never runs with missing/incomplete inputs from an earlier one.
"""

import configparser
import os
import subprocess
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

STEP_ORDER = ("master_bias", "master_dark", "master_flat",
              "calibration", "combination")

COMBINATION_METHODS = ("mean", "median", "sum", "sclip",
                       "medium", "avg", "average")  # incl. aliases

def _master_combine_method(raw):
    """Normalize/validate the master combine method (average/median)."""
    val = {"mean": "average", "avg": "average"}.get(raw.lower(), raw.lower())
    if val not in ("average", "median"):
        raise ValueError(raw)
    return val


# Setting overrides: (section, key) -> (environment variable, caster).
# The step scripts read the CALIB_* variables with their built-in values as
# defaults, so running a script standalone (without the runner) keeps the
# defaults. Keys commented out in config.ini simply don't appear here.
SETTINGS = {
    ("general", "temp_match_tol"):            ("CALIB_TEMP_MATCH_TOL", float),
    ("general", "mem_budget_mb"):             ("CALIB_MEM_BUDGET_MB", int),
    ("master_settings", "combine_method"):    ("CALIB_COMBINE_METHOD",
                                               _master_combine_method),
    ("master_settings", "sigma_clip_low"):    ("CALIB_SIGMA_CLIP_LOW", float),
    ("master_settings", "sigma_clip_high"):   ("CALIB_SIGMA_CLIP_HIGH", float),
    ("master_settings", "min_frames_for_qc"): ("CALIB_MIN_FRAMES_FOR_QC", int),
    ("master_settings", "median_sigma"):      ("CALIB_MEDIAN_SIGMA", float),
    ("master_settings", "median_floor_adu"):  ("CALIB_MEDIAN_FLOOR_ADU", float),
    ("master_settings", "std_ratio_max"):     ("CALIB_STD_RATIO_MAX", float),
    ("master_settings", "max_rn_pairs"):      ("CALIB_MAX_RN_PAIRS", int),
    ("master_dark", "neg_median_tol_adu"):    ("CALIB_NEG_MEDIAN_TOL_ADU", float),
    ("master_flat", "min_frames_for_qc"):     ("CALIB_FLAT_MIN_FRAMES_FOR_QC", int),
    ("master_flat", "saturation_adu"):        ("CALIB_SATURATION_ADU", float),
    ("master_flat", "sat_median_frac"):       ("CALIB_SAT_MEDIAN_FRAC", float),
    ("master_flat", "min_median_frac"):       ("CALIB_MIN_MEDIAN_FRAC", float),
    ("calibration", "flat_min"):              ("CALIB_FLAT_MIN", float),
    ("calibration", "hot_pixel_sigma"):       ("CALIB_HOT_PIXEL_SIGMA", float),
    ("calibration", "cr_sigclip"):            ("CALIB_CR_SIGCLIP", float),
    ("calibration", "cr_objlim"):             ("CALIB_CR_OBJLIM", float),
    ("calibration", "cr_readnoise"):          ("CALIB_CR_READNOISE", float),
}

# Non-numeric per-step options consumed by build_steps, allowed per section.
STEP_OPTIONS = {
    "general": set(),
    "master_settings": set(),
    "master_bias": set(),
    "master_dark": {"no_bias", "level_align"},
    "master_flat": {"dark_scale"},
    "calibration": {"dark_scale", "fix_hot_pixels", "fix_cosmic_rays"},
    "combination": {"method"},
}


def build_env(config):
    """
    Build the subprocess environment from the numeric settings spread over
    the config sections. Every key in a known section must be either a step
    option or a known setting - typos stop the pipeline before it starts.
    Returns (env, applied) - applied is a list of 'section.key=value'
    strings for the startup message.
    """
    if config.has_section("settings") and config.items("settings"):
        print("ERROR: the old [settings] section was replaced - move its "
              "keys into [general] / [master_settings] / [master_dark] / "
              "[master_flat] / [calibration].")
        sys.exit(1)

    env = os.environ.copy()
    applied = []
    for section, step_keys in STEP_OPTIONS.items():
        if not config.has_section(section):
            continue
        for key, raw in config.items(section):
            if key in step_keys:
                continue
            if (section, key) not in SETTINGS:
                valid = sorted(k for s, k in SETTINGS if s == section)
                valid += sorted(step_keys)
                print(f"ERROR: unknown key '{key}' in [{section}]. "
                      f"Available: {', '.join(valid) or '(none)'}")
                sys.exit(1)
            env_name, cast = SETTINGS[(section, key)]
            try:
                value = cast(raw.strip())
            except ValueError:
                print(f"ERROR: [{section}] {key} has an invalid value "
                      f"'{raw}'")
                sys.exit(1)
            env[env_name] = str(value)
            applied.append(f"{section}.{key}={value}")
    return env, applied


def load_config(path: Path):
    if not path.is_file():
        print(f"ERROR: config file not found: {path}")
        sys.exit(1)
    config = configparser.ConfigParser(inline_comment_prefixes=(";", "#"))
    # utf-8-sig: tolerate the BOM that Windows editors often prepend
    config.read(path, encoding="utf-8-sig")
    return config


def build_steps(config):
    """
    Validate the config and return (steps, path_env):
      - steps:    list of (name, command) to run;
      - path_env: CALIB_*_DIR overrides that tell the calibration/combination
                  scripts where to write.

    Layout: everything goes under output_base = <images_path>/<output_folder>
    ("" -> images_path itself). Inside it: <masters_folder>,
    <calibrated_folder>, <combined_folder> (defaults masters/calibrated/
    combined). Any of the four names may be an absolute path, in which case
    it replaces the base (standard pathlib join behaviour).
    """
    try:
        images_path = Path(config.get("paths", "images_path").strip())
    except (configparser.NoSectionError, configparser.NoOptionError) as exc:
        print(f"ERROR: incomplete [paths] section in config: {exc}")
        sys.exit(1)

    output_folder = config.get("paths", "output_folder", fallback="").strip()
    masters_folder = config.get("paths", "masters_folder",
                                fallback="masters").strip() or "masters"
    calibrated_folder = config.get("paths", "calibrated_folder",
                                   fallback="calibrated").strip() or "calibrated"
    combined_folder = config.get("paths", "combined_folder",
                                 fallback="combined").strip() or "combined"

    if not images_path.is_dir():
        print(f"ERROR: images_path does not exist or is not a directory: "
              f"{images_path}")
        sys.exit(1)

    output_base = images_path / output_folder if output_folder else images_path
    masters_dir = output_base / masters_folder
    calibrated_dir = output_base / calibrated_folder
    combined_dir = output_base / combined_folder

    path_env = {
        "CALIB_CALIBRATED_DIR": str(calibrated_dir),
        "CALIB_COMBINED_DIR": str(combined_dir),
    }

    def enabled(step):
        return config.getboolean("steps", step, fallback=False)

    def option(section, key):
        return config.getboolean(section, key, fallback=False)

    steps = []

    if enabled("master_bias"):
        steps.append(("master bias", [
            str(SCRIPT_DIR / "master_bias.py"),
            str(images_path), str(masters_dir),
        ]))

    if enabled("master_dark"):
        cmd = [str(SCRIPT_DIR / "master_dark.py"),
               str(images_path), str(masters_dir)]
        if option("master_dark", "no_bias"):
            cmd.append("no_bias")
        if option("master_dark", "level_align"):
            cmd.append("level_align")
        steps.append(("master dark", cmd))

    if enabled("master_flat"):
        cmd = [str(SCRIPT_DIR / "master_flat.py"),
               str(images_path), str(masters_dir)]
        if option("master_flat", "dark_scale"):
            cmd.append("dark_scale")
        steps.append(("master flat", cmd))

    if enabled("calibration"):
        cmd = [str(SCRIPT_DIR / "light_calibration.py"),
               str(images_path), str(masters_dir)]
        if option("calibration", "dark_scale"):
            cmd.append("dark_scale")
        if option("calibration", "fix_hot_pixels"):
            cmd.append("fix_hot_pixels")
        if option("calibration", "fix_cosmic_rays"):
            cmd.append("fix_cosmic_rays")
        steps.append(("light calibration", cmd))

    if enabled("combination"):
        method = config.get("combination", "method", fallback="mean") \
            .strip().lower() or "mean"
        if method not in COMBINATION_METHODS:
            print(f"ERROR: unknown combination method '{method}' in config "
                  f"(available: {', '.join(COMBINATION_METHODS)})")
            sys.exit(1)
        # light_calibration.py writes into calibrated_dir; when the
        # calibration step runs in this same pipeline the folder will appear,
        # otherwise it must already exist.
        if not enabled("calibration") and not calibrated_dir.is_dir():
            print(f"ERROR: the combination step is enabled but "
                  f"{calibrated_dir} does not exist and the "
                  f"calibration step is disabled. Enable the calibration "
                  f"step or run it first.")
            sys.exit(1)
        steps.append(("light combination", [
            str(SCRIPT_DIR / "light_combination.py"),
            str(calibrated_dir), method,
        ]))

    return steps, path_env


def main():
    config_path = Path(sys.argv[1]) if len(sys.argv) > 1 \
        else SCRIPT_DIR / "config.ini"
    config = load_config(config_path)
    print(f"Config: {config_path}")

    enabled_names = [s for s in STEP_ORDER
                     if config.getboolean("steps", s, fallback=False)]
    if not enabled_names:
        print("No step is enabled in the [steps] section - nothing to do.")
        sys.exit(0)
    print(f"Enabled steps: {', '.join(enabled_names)}")

    steps, path_env = build_steps(config)
    env, applied = build_env(config)
    env.update(path_env)
    if applied:
        print(f"Settings overrides: {', '.join(applied)}")

    t_start = time.time()
    done = []

    for name, cmd in steps:
        print(f"\n{'=' * 70}")
        print(f"STEP: {name}")
        print(f"  {Path(cmd[0]).name} {' '.join(cmd[1:])}")
        print("=" * 70)
        t0 = time.time()
        result = subprocess.run([sys.executable] + cmd, cwd=SCRIPT_DIR,
                                env=env)
        dt = time.time() - t0
        if result.returncode != 0:
            print(f"\n{'=' * 70}")
            print(f"PIPELINE STOPPED: step '{name}' failed "
                  f"(exit code {result.returncode}) after {dt:.1f}s")
            if done:
                print(f"Steps completed before the failure: {', '.join(done)}")
            sys.exit(result.returncode)
        done.append(name)
        print(f"-- step '{name}' finished in {dt:.1f}s")

    print(f"\n{'=' * 70}")
    print(f"PIPELINE DONE: {len(done)} step(s) in {time.time() - t_start:.1f}s "
          f"({', '.join(done)})")


if __name__ == "__main__":
    main()
