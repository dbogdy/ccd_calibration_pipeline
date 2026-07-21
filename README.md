# 📘 CCD Calibration Pipeline

Automated pipeline in Python for **calibrating CCD/CMOS images** using [`astropy`](https://docs.astropy.org/) and [`ccdproc`](https://ccdproc.readthedocs.io/).  
Builds master calibration frames (BIAS, DARK, FLAT) and applies them to science frames (`LIGHT`), then optionally combines the calibrated frames — following the standard CCD reduction procedures from the [Astropy CCD Data Reduction Guide](https://www.astropy.org/ccd-reduction-and-photometry-guide/).

---

## 🧩 Features

- Automatic grouping of FITS files by:
  - `IMAGETYP`
  - binning (`XBINNING`, `YBINNING`)
  - `GAIN`
  - `EXPTIME` (darks / flats / lights)
  - commanded set-point temperature (`SET-TEMP`, exact), with `CCD-TEMP` used as a cooling-stabilization check
  - `FILTER` (flats / lights) and EL injection current (combination)
- Creation of master frames using:
  - sigma-clipped `average` combination (`median` optional), 5σ/5σ by default
  - `mad_std` as the dispersion estimator and `np.ma.median` as the center
  - per-frame quality control (level / noise outliers rejected)
  - read-noise estimate from frame-pair differences (`RDNOISE`)
  - memory-limited, chunked I/O for large stacks
- Science frame calibration (all through `ccdproc`):
  1. Overscan subtraction / trim (when `BIASSEC` / `TRIMSEC` are present)
  2. Bias subtraction (only when the master dark is bias-subtracted)
  3. Dark subtraction (with optional exposure-time scaling)
  4. Flat-field correction with protection against divide-by-zero / dead pixels
  5. Optional hot-pixel repair (from the master dark) and cosmic-ray rejection (L.A.Cosmic)
- Combination of calibrated frames (`mean` / `median` / `sum` / `sclip`)
- Robust handling:
  - automatic master selection based on metadata
  - `no_bias`, `level_align` and `dark_scale` options for tricky CMOS sessions
  - FITS `HISTORY` updated with each processing step
- Fully file-system based, no external database required.

---

## ⚙️ Installation

Requires Python ≥ 3.10.

Core dependencies:
- `numpy`
- `astropy`
- `ccdproc`
- `astroscrappy` (only when cosmic-ray rejection, `fix_cosmic_rays`, is enabled)

Install:

```bash
pip install numpy astropy ccdproc astroscrappy
```

## 📁 Directory Structure
Expected layout before running the pipeline (the `Bias` / `Dark` / `Flat` / `Light` subfolders are matched case-insensitively inside `images_path`):
```bash
images_path/
├── Bias/
├── Dark/
├── Flat/
├── Light/                 # may contain per-object / per-panel subfolders
├── <masters_folder>/      # generated master frames
├── calibrated/            # generated calibrated science frames
└── combined/              # generated combined science frames
```
The pipeline scripts live together (e.g. in the repository root):
```bash
run_pipeline.py            # config-driven runner
config.ini                 # configuration
calib_common.py            # generic helpers (paths, header parsing, formatting)
grouping.py                # frame grouping
masters.py                 # master index + selection
ccd_ops.py                 # ccdproc wrappers, QC and statistics
master_bias.py             # step 1
master_dark.py             # step 2
master_flat.py             # step 3
light_calibration.py       # step 4
light_combination.py       # step 5
```
Input FITS files must:

Contain correct IMAGETYP values: e.g. BIAS, DARK, FLAT, LIGHT

Provide consistent EXPTIME, XBINNING / YBINNING, GAIN

Optionally provide SET-TEMP / CCD-TEMP for temperature-aware grouping, and FILTER for flats/lights

## 🧰 Configuration
The pipeline is controlled via an INI file (default: `config.ini` next to `run_pipeline.py`). It is organized into sections; commented lines equal the built-in defaults, so uncomment only what you want to change.

```ini
[paths]
images_path = /path/to/calibration_files   ; folder with the Bias/Dark/Flat/Light subfolders
masters_folder = processed/masters          ; relative to images_path, or an absolute path

[steps]
master_bias  = true                         ; true / false for each step
master_dark  = true
master_flat  = true
calibration  = true
combination  = true

[general]
temp_match_tol = 1.0                         ; max |CCD-TEMP - SET-TEMP| for a usable frame [deg C]
mem_budget_mb  = 1024                         ; approximate RAM ceiling for combination [MB]

[master_settings]
combine_method    = average                  ; average (recommended) or median
sigma_clip_low    = 5.0                       ; sigma clipping at combination
sigma_clip_high   = 5.0
min_frames_for_qc = 10                        ; below this frame count, QC only warns
median_sigma      = 5.0                       ; reject if |median - group| > median_sigma*MAD + floor
median_floor_adu  = 5.0
std_ratio_max     = 2.0                       ; reject if noise is outside [group/x, group*x]
max_rn_pairs      = 8                         ; frame pairs used for read-noise estimation

[master_dark]
no_bias           = true                      ; build darks from raw frames (recommended for CMOS)
level_align       = true                      ; align each frame's level to the group median
neg_median_tol_adu = 0.2                       ; anti-damage screen against a mismatched master bias

[master_flat]
dark_scale        = false                     ; allow scaling a bias-subtracted dark by exposure
min_frames_for_qc = 5
saturation_adu    = 65535
sat_median_frac   = 0.90
min_median_frac   = 0.05

[calibration]
dark_scale        = false                     ; allow scaling the dark when calibrating lights
flat_min          = 0.2                        ; flat pixels below this are clamped (min_value)
fix_cosmic_rays   = false                     ; L.A.Cosmic rejection per frame (needs astroscrappy)
cr_sigclip        = 4.5                        ; L.A.Cosmic detection threshold [sigma]
cr_objlim         = 5.0                        ; L.A.Cosmic CR-vs-source contrast limit
cr_readnoise      = 8.0                        ; fallback read noise [e-] when RDNOISE is absent
fix_hot_pixels    = false                     ; interpolate hot pixels flagged from the master dark
hot_pixel_sigma   = 5.0                        ; hot-pixel threshold [sigma above the dark median]

[combination]
method            = mean                      ; mean / median / sum / sclip
```

`run_pipeline.py` passes the numeric settings to the step scripts through `CALIB_*` environment variables; a script run standalone (without the runner) uses the built-in defaults from the code. An unknown key in a known section stops the pipeline with an error (typo protection).

## ▶️ Usage

Run the whole pipeline using the default configuration file:
```bash
    python run_pipeline.py
```
or specify a custom one:
```bash
    python run_pipeline.py path/to/config.ini
```
Each step can also be run standalone:
```bash
    python master_bias.py        <images_path> [masters_folder]
    python master_dark.py        <images_path> [masters_folder] [no_bias] [level_align]
    python master_flat.py        <images_path> [masters_folder] [dark_scale]
    python light_calibration.py  <images_path> [masters_folder] [dark_scale]
    python light_combination.py  <images_path>/calibrated [method]
```

## 📊 Output

| Type                      | Directory      | Example Filename                                        |
| ------------------------- | -------------- | ------------------------------------------------------ |
| Master Bias               | `<masters>/`   | `master_bias_bin2x2_gain120_temp-20C.fits`             |
| Master Dark               | `<masters>/`   | `master_dark_bin2x2_gain120_exp10s_temp-20C.fits`      |
| Master Flat               | `<masters>/`   | `master_flat_bin2x2_gain120_exp1s_temp-20C_filtSDSS-y.fits` |
| Calibrated Science Frames | `calibrated/`  | `Light_001_cal.fits`                                   |
| Combined Science Frames   | `combined/`    | `combined_mean_bin2x2_gain120_exp10s_temp-20C_5A_SDSS-y.fits` |

## How it works (technical summary)

1. Master builders (`master_bias.py`, `master_dark.py`, `master_flat.py`)

    - Scan the `Bias` / `Dark` / `Flat` subfolders for FITS files.

    - Group frames by `(IMAGETYP, binning, GAIN, EXPTIME, SET-TEMP)`; frames whose `CCD-TEMP` deviates from `SET-TEMP` by more than `temp_match_tol` are dropped (cooler not stabilized).

    - Run per-frame QC, then build master BIAS/DARK/FLAT with `ccdproc.combine` (sigma clipping, `mad_std`); flats are combined with `scale=1/median` and renormalized to a median of 1.0.

    - Encode the grouping parameters in the output filename and update the FITS `HISTORY`.

2. Science calibration (`light_calibration.py`)

    - Groups the `Light` frames, selects the matching masters and calibrates each frame with `ccdproc`: overscan/trim → `subtract_bias` (only when the dark is bias-subtracted) → `subtract_dark` (scaled when allowed) → `flat_correct`.

    - Optionally cleans each calibrated frame: hot-pixel interpolation (from the master dark) and L.A.Cosmic cosmic-ray rejection, both off by default and enabled via `[calibration]`.

    - Writes `*_cal.fits` into `calibrated/`, preserving the Light subfolder structure.

3. Combination (`light_combination.py`)

    - Groups the calibrated frames (subfolder, binning, gain, exposure, temperature, `FILTER`, EL current) and combines each multi-frame group with the chosen method into `combined/`.

    - `FILTER` and EL current are read from the header, or parsed from the filename when absent.

4. Runner (`run_pipeline.py`)

    - Reads `config.ini`, validates it, and runs the enabled steps in order, stopping at the first failing step so later steps never run with missing inputs.

## Assumptions and limitations

- Assumes primary image data in the primary HDU.

- Assumes standard FITS keywords:

    - IMAGETYP, EXPTIME, XBINNING, YBINNING, GAIN, SET-TEMP / CCD-TEMP (and FILTER for flats/lights).

- Designed for single-CCD workflows; multi-extension instruments need adaptation.

- Relies on correct header metadata for reliable master selection.

## 🔬 Scientific Notes

- Master combination uses sigma clipping (5σ by default, configurable) with `mad_std` dispersion around the median; `sclip` in the combination step uses 3σ.

- Temperature matching is exact on `SET-TEMP` (the commanded set point is discrete, so no tolerance chaining is needed); `CCD-TEMP` is only a stabilization check with tolerance `temp_match_tol` (default 1 °C). Frames without `SET-TEMP` fall back to the rounded `CCD-TEMP`.

- Exposure matching uses a relative tolerance of 1e-6.

- A dark is scaled by the exposure ratio only when it is bias-subtracted; a median-level screen blocks a master bias that would drive the darks negative (CMOS with an unstable bias level).

- Optional per-frame cleaning (both off by default): `fix_hot_pixels` flags hot pixels from the master dark (`clipped_median + hot_pixel_sigma·clipped_std`) and replaces them with the 3×3 neighbour median; `fix_cosmic_rays` runs L.A.Cosmic. These are cosmetic per-frame steps — for pure statistical outlier rejection across a stack, prefer the `sclip` combination method.

- FITS metadata and processing history are preserved in the output headers.

## 🧠 Future Improvements

- Optional bad-pixel map to mask dead / over-corrected pixels

- Multi-extension (multi-CCD) support

## Acknowledgement
- This work was supported by the Romanian Ministry of Education and Research, through the Executive Agency for Higher Education, Research, Development and Innovation Funding (UEFISCDI), under the National Plan for Research, Development and Innovation 2022–2027 (PN IV), Demonstration Project (PED), Contract No. 12PED/2025, project "IMAGINER – Image Enhancement Algorithms for Photovoltaic Panel Monitoring".

## References

- Astropy Project: core library and documentation (https://www.astropy.org/)

- Astropy CCD Data Reduction Guide (https://www.astropy.org/ccd-reduction-and-photometry-guide/)

- ccdproc: CCD data reduction tools (https://ccdproc.readthedocs.io/)
                                    (https://github.com/astropy/ccdproc)
