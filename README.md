# 📘 CCD Calibration Pipeline

A Python-based automated pipeline for **calibrating astronomical CCD images**, using `astropy` and `ccdproc`.  
It builds master calibration frames (BIAS, DARK, FLAT) and applies them to science frames (LIGHT) following standard CCD data reduction practices.
---
## 🧩 Features

- Automatic grouping of FITS files by metadata (`IMAGETYP`, `XBINNING`, `GAIN`, `EXPTIME`, `CCD-TEMP`)
- Master frame creation using sigma-clipped median combining with robust dispersion (`mad_std`)
- Full calibration workflow:
  1. Bias subtraction  
  2. Dark subtraction (optional scaling by exposure time)  
  3. Flat-field correction (adaptive normalization)
- Robust numerical handling (no divide-by-zero errors)
- Cross-matching of masters by exposure and temperature tolerance
- FITS headers store full `HISTORY` of each processing step
---
## ⚙️ Requirements

| Dependency  | Version  | Purpose                    |
|-------------|----------|----------------------------|
| Python      |  ≥3.10   | Core runtime               |
| numpy       |  latest  | array operations           |
| astropy     |   ≥6.0   | FITS I/O and units         |
| ccdproc     |  latest  | image calibration routines |

Install manually:
```bash
    pip install numpy astropy ccdproc
```
or use 
```bash
    pip install -r requirements.txt
```

📁 Directory Structure
Expected layout before running the pipeline:

project_root/
│
├── raw/
│   ├── BIAS/
│   ├── DARK/
│   ├── FLAT/
│   └── LIGHT/
│
└──calib_dir
    ├── masters/        # generated master frames
    ├── calibrated/     # output calibrated images
    ├── config.ini
    ├── fits_image_calibration.py
    ├── run_calibration_pipeline.py
    └── README.md

🧰 Configuration (config.ini)
Example:
[GENERAL]
debug = yes

[MASTERS]
raw_root = ../raw
master_output = masters
build_bias = yes
build_dark = yes
build_flat = yes
dark_scale = no

[CALIBRATION]
input_root = ../raw
output_dir = calibrated
master_files = masters
file_type = LIGHT
no_bias = no
dark_scale = no

| Section         | Option       | Description                           |
| --------------- | ------------ | ------------------------------------- |
| `[GENERAL]`     | `debug`      | Enables verbose console logging       |
| `[MASTERS]`     | `build_*`    | Control which master frames are built |
| `[CALIBRATION]` | `output_dir` | Output folder for calibrated FITS     |

▶️ Usage

Run the pipeline using the configuration file:
```bash
    python run_calibration_pipeline.py
```
or specify a custom one:
```bash
    python run_calibration_pipeline.py path/to/config.ini
```
📊 Output

| Type                      | Directory     | Example Filename                      |
| ------------------------- | ------------- | ------------------------------------- |
| Master Bias               | `masters/`    | `master_bias_b2x2_g120_sNA_C-20.fits` |
| Master Dark               | `masters/`    | `master_dark_b2x2_g120_10s_C-20.fits` |
| Master Flat               | `masters/`    | `master_flat_b2x2_g120_sNA_C-20.fits` |
| Calibrated Science Frames | `calibrated/` | `Light_001_cal.fits`                  |

🔬 Scientific Notes

Combination uses SigmaClip (3σ, 5 iterations) with mad_std dispersion

Master matching tolerates ±2 °C in temperature and relative exposure difference ≤ 1e-6

Bias and dark corrections optional and fully configurable

FITS metadata and processing history preserved in output headers

🧠 Future Improvements

Automatic filter keyword handling (FILTER/FILT)
