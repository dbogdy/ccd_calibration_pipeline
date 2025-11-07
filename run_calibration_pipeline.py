#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Runner for CCD master creation and science frame calibration.

Configuration is read from an INI file (default: config.ini).
See example in the repository root.

Usage
-----
    python run_ccd_pipeline.py            # uses ./config.ini
    python run_ccd_pipeline.py my.cfg     # uses custom config file
"""

import sys
import configparser
from pathlib import Path

from fits_image_calibration import create_master_file, calibrate_files


def _get_bool(cfg: configparser.ConfigParser, section: str, option: str, default: bool = False) -> bool:
    """
    Read a boolean-like option from config.

    Accepts: yes/no, true/false, 1/0. Case-insensitive.
    """
    if not cfg.has_option(section, option):
        return default
    return cfg.get(section, option).strip().lower() in {"1", "true", "yes", "y"}


def main(config_path: str = "config.ini") -> None:
    """
    Main entry point.

    Steps
    -----
    1. Read configuration.
    2. Optionally build BIAS, DARK, FLAT masters.
    3. Optionally calibrate science frames using the created masters.
    """
    cfg = configparser.ConfigParser()
    cfg.read(config_path)

    # General options
    debug = _get_bool(cfg, "GENERAL", "debug", default=False)

    # Master frame options
    raw_root = Path(cfg.get("MASTERS", "raw_root", fallback="./raw")).resolve()
    master_output = Path(cfg.get("MASTERS", "master_output", fallback="./masters")).resolve()

    build_bias = _get_bool(cfg, "MASTERS", "build_bias", default=True)
    build_dark = _get_bool(cfg, "MASTERS", "build_dark", default=True)
    build_flat = _get_bool(cfg, "MASTERS", "build_flat", default=True)

    dark_no_bias = _get_bool(cfg, "MASTERS", "dark_no_bias", default=False)
    dark_scale = _get_bool(cfg, "MASTERS", "dark_scale", default=False)

    flat_no_bias = _get_bool(cfg, "MASTERS", "flat_no_bias", default=False)
    flat_dark_scale = _get_bool(cfg, "MASTERS", "flat_dark_scale", default=False)

    # Calibration options
    file_type = cfg.get("CALIBRATION", "file_type", fallback="LIGHT").strip().upper()
    input_root = Path(cfg.get("CALIBRATION", "input_root", fallback="./raw")).resolve()
    output_dir = Path(cfg.get("CALIBRATION", "output_dir", fallback="./calibrated")).resolve()

    master_files_str = cfg.get("CALIBRATION", "master_files", fallback=str(master_output))
    master_files = [p.strip() for p in master_files_str.split(",") if p.strip()]

    cal_no_bias = _get_bool(cfg, "CALIBRATION", "no_bias", default=False)
    cal_dark_scale = _get_bool(cfg, "CALIBRATION", "dark_scale", default=False)

    # 1) Build master BIAS
    if build_bias:
        create_master_file(
            dir_path=raw_root,
            file_type="BIAS",
            output=master_output,
            master_files=None,
            no_bias=True,          # Bias should not be bias-corrected
            dark_scale=False,
            debug=debug,
            mem_limit=None,
        )

    # 2) Build master DARK
    if build_dark:
        create_master_file(
            dir_path=raw_root,
            file_type="DARK",
            output=master_output,
            master_files=str(master_output),
            no_bias=dark_no_bias,
            dark_scale=dark_scale,
            debug=debug,
            mem_limit=None,
        )

    # 3) Build master FLAT
    if build_flat:
        create_master_file(
            dir_path=raw_root,
            file_type="FLAT",
            output=master_output,
            master_files=str(master_output),
            no_bias=flat_no_bias,
            dark_scale=flat_dark_scale,
            debug=debug,
            mem_limit=None,
        )

    # 4) Calibrate science frames
    if cfg.has_section("CALIBRATION"):
        calibrate_files(
            input_dir=input_root,
            file_type=file_type,
            output_dir=output_dir,
            master_files=master_files,
            no_bias=cal_no_bias,
            dark_scale=cal_dark_scale,
            debug=debug,
        )


if __name__ == "__main__":
    # Allow optional custom config path as first CLI argument.
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main()
