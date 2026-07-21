#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared low-level helpers for the ccdproc-based calibration pipeline (v3).

Only generic, non-pixel logic lives here: path lookup, FITS header value
parsing and small formatting helpers. All pixel processing is done with
ccdproc (see ccd_ops.py).
"""

import math
import os
from pathlib import Path


def env_float(name, default):
    """
    Read a setting override from the environment. run_pipeline.py sets these
    variables from the config.ini sections ([general] / [master_settings] /
    [master_dark] / [master_flat] / [calibration]); scripts run standalone
    (no runner) keep the built-in default.
    """
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return float(raw)
    except ValueError:
        print(f"WARNING: ignoring invalid value '{raw}' for {name}, "
              f"using default {default}")
        return default


def env_int(name, default):
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return int(raw)
    except ValueError:
        print(f"WARNING: ignoring invalid value '{raw}' for {name}, "
              f"using default {default}")
        return default


def env_choice(name, default, choices):
    raw = os.environ.get(name)
    if raw is None:
        return default
    val = raw.strip().lower()
    if val in choices:
        return val
    print(f"WARNING: ignoring invalid value '{raw}' for {name} "
          f"(allowed: {', '.join(choices)}), using default {default}")
    return default


# Resources
MEM_BUDGET_MB = env_int("CALIB_MEM_BUDGET_MB", 1024)
# approx. RAM ceiling for ccdproc.combine


def find_subfolder(root: Path, name: str):
    """Return the subfolder of root with this name (case-insensitive), or None."""
    for entry in root.iterdir():
        if entry.is_dir() and entry.name.lower() == name.lower():
            return entry
    return None


def collect_fits_files(folder: Path):
    """All .fit/.fits files in folder (recursive, case-insensitive extension)."""
    return sorted(
        f for f in folder.rglob("*")
        if f.is_file() and f.suffix.lower() in (".fit", ".fits")
    )


def as_int(value):
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return None


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def exposures_match(a, b):
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    return math.isclose(a, b, rel_tol=1e-6)


def fexp(v):
    """Format an exposure value for messages and file names ('NA' when None)."""
    return f"{v:g}" if v is not None else "NA"


def tag(value, prefix=""):
    """Format a header value for file names ('NA' when None)."""
    return f"{prefix}{value}" if value is not None else f"{prefix}NA"
