#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Runner for CCD master creation and science frame calibration.

Configuration is read from an INI file (default: config.ini).
See example in the repository root.

Usage
-----
    python run_calibration_pipeline.py            # uses ./config.ini
    python run_calibration_pipeline.py my.cfg     # uses custom config file
"""

import logging
import sys
import configparser
from pathlib import Path
from typing import Union

from fits_image_calibration import create_master_file, calibrate_files

def _get_bool(cfg, section, option, default=False):
    try:
        return cfg.getboolean(section, option)
    except (configparser.NoSectionError, configparser.NoOptionError, ValueError):
        return default

def _resolve_cfg_path(cfg, section, option, fallback, base_dir: Path) -> Path:
    """
    Reads a path from config and, if it is relative,
    makes it relative to the config file directory.
    """
    val = cfg.get(section, option, fallback=fallback).strip()
    p = Path(val)
    if not p.is_absolute():
        p = (base_dir / p).resolve()
    else:
        p = p.resolve()
    return p

def main(config_path: Union[str, Path] = "config.ini") -> None:
    """
    Main entry point.

    Steps
    -----
    1. Read configuration.
    2. Optionally build BIAS, DARK, FLAT masters.
    3. Optionally calibrate science frames using the created masters.
    """

    # ---- Find the config file directory ----
    config_path = Path(config_path).resolve()
    config_dir = config_path.parent

    cfg = configparser.ConfigParser()
    cfg.read(str(config_path))

    # General options
    debug = _get_bool(cfg, "GENERAL", "debug", default=False)

    # ---------- Master frame options ----------
    raw_root = _resolve_cfg_path(cfg, "MASTERS", "raw_root", "./raw", config_dir)
    master_output = _resolve_cfg_path(cfg, "MASTERS", "master_output", "./masters", config_dir)

    build_bias = _get_bool(cfg, "MASTERS", "build_bias", default=True)
    build_dark = _get_bool(cfg, "MASTERS", "build_dark", default=True)
    build_flat = _get_bool(cfg, "MASTERS", "build_flat", default=True)

    dark_no_bias = _get_bool(cfg, "MASTERS", "dark_no_bias", default=False)
    dark_scale   = _get_bool(cfg, "MASTERS", "dark_scale", default=False)

    flat_no_bias    = _get_bool(cfg, "MASTERS", "flat_no_bias", default=False)
    flat_dark_scale = _get_bool(cfg, "MASTERS", "flat_dark_scale", default=False)

    mem_limit_mb = cfg.getint("MASTERS", "mem_limit_mb", fallback=2048)
    mem_limit = mem_limit_mb * 1024 * 1024

    # ---------- Calibration options ----------
    run_calibration = _get_bool(cfg, "CALIBRATION", "run_calibration", default=True)

    file_type = cfg.get("CALIBRATION", "file_type", fallback="LIGHT").strip().upper()
    input_root = _resolve_cfg_path(cfg, "CALIBRATION", "input_root", "./raw", config_dir)
    output_dir = _resolve_cfg_path(cfg, "CALIBRATION", "output_dir", "./calibrated", config_dir)

    master_files_str = cfg.get("CALIBRATION", "master_files", fallback=str(master_output))
    master_files_raw = [p.strip() for p in master_files_str.split(",") if p.strip()]

    # resolve master_files paths relative to config_dir
    master_files = []
    for s in master_files_raw:
        p = Path(s)
        if not p.is_absolute():
            p = (config_dir / p).resolve()
        else:
            p = p.resolve()
        master_files.append(str(p))

    cal_no_bias    = _get_bool(cfg, "CALIBRATION", "no_bias", default=False)
    cal_dark_scale = _get_bool(cfg, "CALIBRATION", "dark_scale", default=False)

    fix_hot_pixels  = _get_bool(cfg, "CALIBRATION", "fix_hot_pixels",  default=True)
    hot_pixel_sigma = cfg.getfloat("CALIBRATION", "hot_pixel_sigma", fallback=5.0)

    fix_cosmic_rays = _get_bool(cfg, "CALIBRATION", "fix_cosmic_rays", default=True)
    cr_sigclip      = cfg.getfloat("CALIBRATION", "cr_sigclip",      fallback=4.5)
    cr_objlim       = cfg.getfloat("CALIBRATION", "cr_objlim",       fallback=5.0)
    cr_readnoise    = cfg.getfloat("CALIBRATION", "cr_readnoise",    fallback=8.0)

    # 1) Build master BIAS
    if build_bias:
        create_master_file(
            dir_path=raw_root,
            file_type="BIAS",
            output=master_output,
            master_files=None,
            no_bias=True,
            dark_scale=False,
            debug=debug,
            mem_limit=mem_limit,
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
            mem_limit=mem_limit,
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
            mem_limit=mem_limit,
        )

    # 4) Calibrate science frames
    if cfg.has_section("CALIBRATION") and run_calibration:
        calibrate_files(
            input_dir=input_root,
            file_type=file_type,
            output_dir=output_dir,
            master_files=master_files,
            no_bias=cal_no_bias,
            dark_scale=cal_dark_scale,
            debug=debug,
            fix_hot_pixels=fix_hot_pixels,
            hot_pixel_sigma=hot_pixel_sigma,
            fix_cosmic_rays=fix_cosmic_rays,
            cr_sigclip=cr_sigclip,
            cr_objlim=cr_objlim,
            cr_readnoise=cr_readnoise,
        )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
        stream=sys.stderr,
    )
    # Allow optional custom config path as first CLI argument.
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main()
