# CCD Calibration Pipeline

Automated pipeline in Python for **calibrating CCD images** using [`astropy`](https://docs.astropy.org/) and [`ccdproc`](https://ccdproc.readthedocs.io/).
Builds master calibration frames (BIAS, DARK, FLAT) and applies them to science frames (e.g. `LIGHT`) following standard CCD reduction procedures.

---

## Features

- Automatic grouping of FITS files by header keywords:
  - `IMAGETYP`
  - binning (`XBINNING`, `YBINNING`)
  - `GAIN`
  - `EXPTIME`
  - rounded `CCD-TEMP` with ±2 °C tolerance-based matching
  - `FILTER` (for FLAT and LIGHT frames; defaults to `"Clear"` when absent)
- All metadata read exclusively from FITS headers — no filename parsing
- Master frame creation:
  - sigma-clipped median combination (`mad_std` dispersion estimator)
  - optional bias and dark pre-correction for DARK/FLAT masters
- Science frame calibration sequence per group:
  1. Bias subtraction (optional)
  2. Dark subtraction matched by exposure time and temperature (with optional scaling)
  3. Flat-field correction with dead-pixel floor protection
  4. Hot pixel interpolation from master dark (sigma-clipped detection, 3×3 neighbour median)
  5. Cosmic ray rejection via L.A.Cosmic algorithm (optional; uses `EGAIN` header for noise model)
- Robust master selection:
  - ambiguous matches resolved interactively via stdin
  - FITS `HISTORY` updated at each processing step
  - configurable logging with `debug` mode
- Output organised as `calibrated/<object>/<filter>/`
- Fully file-system based, no external database required

---

## Installation

Requires Python ≥ 3.10.

Dependencies:
```
numpy
astropy
ccdproc
```

```bash
pip install numpy astropy ccdproc
```

---

## Directory Structure

Expected layout before running:

```
project_root/
├── raw/
│   ├── BIAS/
│   ├── DARK/
│   ├── FLAT/
│   └── LIGHT/
│       └── ObjectName/       # subdirectory name becomes the object label
├── masters/                  # generated master frames (configured via master_output)
├── calibrated/               # calibrated science frames (configured via output_dir)
│   └── ObjectName/
│       └── FilterName/
├── config.ini
├── fits_image_calibration.py
├── run_calibration_pipeline.py
└── README.md
```

Input FITS files must:
- Contain correct `IMAGETYP` values: `BIAS`, `DARK`, `FLAT`, `LIGHT`
- Provide consistent `EXPTIME`, `XBINNING` / `YBINNING`, `GAIN`
- Optionally provide `CCD-TEMP` for temperature-aware grouping and master matching
- Optionally provide `FILTER` for flat/light grouping (defaults to `"Clear"` if absent)
- Optionally provide `EGAIN` [e⁻/ADU] for accurate cosmic ray noise modelling

---

## Configuration

The pipeline is controlled via an INI file (default: `config.ini`).
All paths can be absolute or relative to the config file location.

```ini
[GENERAL]
debug = yes        # enable verbose DEBUG logging
```

```ini
[MASTERS]
raw_root = ./raw           # root folder containing BIAS/DARK/FLAT subfolders
master_output = ./masters  # where master files are written

build_bias = yes
build_dark = yes
build_flat = yes

dark_no_bias = no          # yes: skip bias subtraction when building master DARK
dark_scale = no            # yes: scale darks by exposure time

flat_no_bias = no          # yes: skip bias subtraction when building master FLAT
flat_dark_scale = no       # yes: scale darks when correcting flats

mem_limit_mb = 2048        # memory limit for frame combination (MB)
```

```ini
[CALIBRATION]
run_calibration = yes
file_type = LIGHT          # IMAGETYP of science frames to calibrate
input_root = ./raw         # root directory containing <file_type> subfolder
output_dir = ./calibrated  # root output directory

# directory or comma-separated list of master FITS paths
master_files = ./masters

no_bias = no               # yes: skip bias subtraction on science frames
dark_scale = no            # yes: scale master darks by exposure time

# Hot pixel correction
# Pixels in the master dark above  median + hot_pixel_sigma * clipped_std
# are flagged and replaced by the median of valid 3×3 neighbours in each LIGHT frame.
fix_hot_pixels = yes
hot_pixel_sigma = 5.0

# Cosmic ray rejection (L.A.Cosmic algorithm)
# Applied after bias+dark+flat+hot-pixel correction.
# cr_readnoise is used only when RDNOISE is absent from the FITS header.
# NOTE: not suitable for images of extended bright objects (e.g. illuminated panels);
#       set fix_cosmic_rays = no in such cases.
fix_cosmic_rays = yes
cr_sigclip   = 4.5
cr_objlim    = 5.0
cr_readnoise = 8.0
```

---

## Usage

```bash
python run_calibration_pipeline.py                  # uses ./config.ini
python run_calibration_pipeline.py path/to/my.cfg   # custom config file
```

---

## Output

| Type | Directory | Example filename |
|------|-----------|-----------------|
| Master Bias | `masters/` | `master_bias_b2x2_g120_sNA_CNA.fits` |
| Master Dark | `masters/` | `master_dark_b2x2_g120_10s_-15.0C.fits` |
| Master Flat | `masters/` | `master_flat_b2x2_g120_1s_-15.0C_Clear.fits` |
| Calibrated science | `calibrated/<object>/<filter>/` | `Light_001_cal.fits` |

---

## How It Works

### `fits_image_calibration.py`

1. Scans input directories recursively for `.fit` / `.fits` files (case-insensitive).
2. Groups frames by `(IMAGETYP, binning, GAIN, EXPTIME, rounded CCD-TEMP, FILTER)`.
   - Temperature rounding step: 0.5 °C; merge tolerance: ±2 °C.
   - FILTER applied to FLAT and LIGHT only; BIAS/DARK always get `None`.
3. Builds master BIAS/DARK using `ccdproc.combine` (sigma-clipped median, `mad_std`).
4. Builds master FLAT using per-frame 1/median scaling before combination.
5. For science frames, finds the best master for each group:
   - DARK: must match exposure exactly (rel_tol=1e-6) and temperature (±2 °C).
   - FLAT: must match filter; exposure and temperature not required.
   - BIAS: matched by binning and gain only.
   - Ambiguous matches prompt the user interactively.
6. Applies calibration: bias → dark → flat → hot pixels → cosmic rays.
7. Writes output as float32 FITS with updated `HISTORY` keywords.

### `run_calibration_pipeline.py`

- Reads configuration from `config.ini`.
- Optionally builds BIAS, DARK, FLAT masters in order.
- Calibrates science frames using the produced (or pre-existing) masters.

---

## Assumptions and Limitations

- Primary image data must be in the primary HDU (extension 0).
- Standard FITS keywords required: `IMAGETYP`, `EXPTIME`, `XBINNING`, `YBINNING`, `GAIN`.
- For accurate cosmic ray noise modelling, `EGAIN` [e⁻/ADU] should be present in the header (falls back to `GAIN` if absent).
- Designed for single-CCD workflows; multi-extension instruments need adaptation.
- L.A.Cosmic is designed for astronomical images with sparse point sources and dark backgrounds. It is not suitable for images of bright extended objects (solar panels, lab targets, etc.) — disable with `fix_cosmic_rays = no` in those cases.

---

## Scientific Notes

- Combination uses sigma-clipping (3σ, 5 iterations) with `mad_std` dispersion.
- Hot pixel detection uses `sigma_clipped_stats` (sigma=3.0) so the threshold is robust even when the pixel distribution is highly quantized (MAD = 0).
- Dark master matching requires exact exposure time (rel_tol=1e-6) and temperature within ±2 °C.
- `dark_scale = yes` linearly scales the dark by the exposure ratio — intended only for cases where an exactly-matched dark is unavailable.
- FITS metadata and processing history preserved in all output headers.

---

## References

- Astropy Project: [https://www.astropy.org/](https://www.astropy.org/)
- ccdproc documentation: [https://ccdproc.readthedocs.io/](https://ccdproc.readthedocs.io/)
- van Dokkum (2001) — L.A.Cosmic algorithm: [https://www.astro.yale.edu/dokkum/lacosmic/](https://www.astro.yale.edu/dokkum/lacosmic/)
