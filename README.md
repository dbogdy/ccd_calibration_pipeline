# 📘 CCD Calibration Pipeline

Automated pipeline in Python for **calibrating CCD images** using [`astropy`](https://docs.astropy.org/) and [`ccdproc`](https://ccdproc.readthedocs.io/).  
Builds master calibration frames (BIAS, DARK, FLAT) and applies them to science frames (e.g. `LIGHT`) following standard CCD reduction procedures.

---

## 🧩 Features

- Automatic grouping of FITS files by:
  - `IMAGETYP`
  - binning (`XBINNING`, `YBINNING`)
  - `GAIN`
  - `EXPTIME`
  - rounded `CCD-TEMP` with tolerance-based matching
- Creation of master frames using:
  - sigma-clipped median combination
  - `mad_std` as dispersion estimator
  - optional bias and dark pre-correction for DARK/FLAT/LIGHT
- Science frame calibration:
  1. Bias subtraction (optional)
  2. Dark subtraction (with optional exposure-time scaling)
  3. Flat-field correction with protection against divide-by-zero / dead pixels
- Robust handling:
  - automatic master selection based on metadata
  - FITS `HISTORY` updated with each processing step
  - logging support (`debug` mode)
- Fully file-system based, no external database required.

---

## ⚙️ Installation

Requires Python ≥ 3.10.

Core dependencies:
- `numpy`
- `astropy`
- `ccdproc`

Install:

```bash
pip install numpy astropy ccdproc
```


## 📁 Directory Structure
Expected layout before running the pipeline:
```bash
project_root/
├── raw/
│   ├── BIAS/
│   ├── DARK/
│   ├── FLAT/
│   └── LIGHT/              # or other science type matching IMAGETYP
├── calibration 
│   ├── masters/            # generated master frames
│   └── calibrated/         # calibrated science frames
├── config.ini
├── fits_image_calibration.py
├── run_calibration_pipeline.py
└── README.md
```
Input FITS files must:

Contain correct IMAGETYP values: e.g. BIAS, DARK, FLAT, LIGHT

Provide consistent EXPTIME, XBINNING / YBINNING, GAIN

Optionally provide CCD-TEMP for temperature-aware grouping

## 🧰 Configuration
The pipeline is controlled via an INI file (default: config.ini).

Example:
```ini
[GENERAL]
debug = yes        # enable verbose logging for internal operations
```

```ini
[MASTERS]
raw_root = ./raw          # root folder containing BIAS/DARK/FLAT subfolders
master_output = ./masters # where master files are written

build_bias = yes
build_dark = yes
build_flat = yes

dark_no_bias = no         # if yes: do not bias-correct darks when building master DARK
dark_scale = no           # if yes: scale darks by exposure time

flat_no_bias = no         # if yes: do not bias-correct flats
flat_dark_scale = no      # if yes: scale darks when correcting flats
```

```ini
[CALIBRATION]
file_type = LIGHT         # IMAGETYP for science frames to calibrate
input_root = ./raw        # root directory containing science subfolder (e.g. ./raw/LIGHT)
output_dir = ./calibrated # where calibrated FITS files will be stored

# can be a directory or comma-separated list of FITS master paths
master_files = ./masters

no_bias = no              # if yes: skip bias subtraction on science frames
dark_scale = no           # if yes: scale master darks by exposure time
```

## ▶️ Usage

Run the pipeline using the configuration file:
```bash
    python run_calibration_pipeline.py
```
or specify a custom one:
```bash
    python run_calibration_pipeline.py path/to/config.ini
```

## 📊 Output

| Type                      | Directory     | Example Filename                      |
| ------------------------- | ------------- | ------------------------------------- |
| Master Bias               | `masters/`    | `master_bias_b2x2_g120_sNA_C-20.fits` |
| Master Dark               | `masters/`    | `master_dark_b2x2_g120_10s_C-20.fits` |
| Master Flat               | `masters/`    | `master_flat_b2x2_g120_sNA_C-20.fits` |
| Calibrated Science Frames | `calibrated/` | `Light_001_cal.fits`                  |

## How it works (technical summary)

1. fits_image_calibration.py

    - Scans input directories for FITS files.

    - Groups frames by (IMAGETYP, binning, GAIN, EXPTIME, rounded CCD-TEMP).

    - Builds master BIAS/DARK/FLAT using ccdproc.combine with sigma clipping and mad_std.

    - Encodes grouping parameters in output filenames.

    - Preserves and updates FITS HISTORY.

2. run_calibration_pipeline.py

    - Cite configuration from config.ini.

    - Optionally builds master frames.

    - Calibrates requested science frames with selected masters.


## Assumptions and limitations

- Assumes primary image data in the primary HDU.

- Assumes standard FITS keywords:

    - IMAGETYP, EXPTIME, XBINNING, YBINNING, GAIN, CCD-TEMP.

- Designed for single-CCD workflows; multi-extension instruments need adaptation.

- Relies on correct header metadata for reliable master selection.

## 🔬 Scientific Notes

- Combination uses SigmaClip (3σ, 5 iterations) with mad_std dispersion

- Master matching tolerates ±2 °C in temperature and relative exposure difference ≤ 1e-6

- Bias and dark corrections optional and fully configurable

- FITS metadata and processing history preserved in output headers

## 🧠 Future Improvements

- Automatic filter keyword handling (FILTER/FILT)

## References

- Astropy Project: core library and documentation (https://www.astropy.org/)

- Astropy documentation (https://docs.astropy.org/)

- ccdproc: CCD data reduction tools (https://ccdproc.readthedocs.io/)
                                    (https://github.com/astropy/ccdproc)


