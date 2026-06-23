# detonation-structure-extraction

Detonation reaction-zone structure analysis from schlieren / shadowgraph image
sequences. Given a sequence of frames of a detonation wave, the tool detects the
wavefront in each frame, aligns every frame to that front, and builds a
coverage-corrected map of where structural features sit relative to the front —
plus the reaction-zone width distribution and cross-dataset comparisons.

Python port of `reactionZone_v5.m` / `reactionZone_v6.m`. Deliberate behavior
changes from the MATLAB original are recorded in `CHANGELOG_PORT.md`.

## What it does

For one dataset (a folder of frames, with frame `001` as the background):

1. Subtract the background frame from every frame and zero-median normalize.
2. Segment structural pixels with a folded-Otsu (or noise-floor) threshold.
3. Detect the wavefront as the rightmost segmented pixel per row, then smooth it.
4. Shift every row so the front lands near the right edge — this maximizes the
   burnt-gas region that stays in view.
5. Track per-pixel *coverage* (which pixels hold real data vs. alignment padding)
   so the final probabilities are normalized by real data, not raw frame count.
6. Filter out frames whose front sits at the field-of-view boundary.
7. Accumulate the aligned masks in randomized order into a probability heatmap,
   and track how the probability at fixed distances converges as frames are added.
8. Build the reaction-zone width distribution (front → last structural pixel per
   row), with censoring of rows truncated by the field of view.

## Requirements

- Python ≥ 3.12
- [`uv`](https://docs.astral.sh/uv/) for environment / dependency management
- Core dependencies (installed automatically by `uv`): `numpy`, `scipy`,
  `scikit-image`, `matplotlib`, `imageio`, and `cantera`
- **Optional** chemistry layer (`chemistry.py`) additionally needs Cantera *and*
  SDToolbox on your path. You only need this if you want CJ velocity / induction
  length computed from a reaction mechanism instead of typing measured values.

## Install

```bash
uv sync          # creates .venv and installs from pyproject.toml / uv.lock
```

Run scripts with `uv run` so the project environment is used:

```bash
uv run python run_analysis.py
```

## Project layout

```
detrz_extraction/        package (importable modules)
  parameters.py            ExtractionParams + PhysicalParams (all the knobs)
  background.py            frame discovery, loading, background subtraction
  segmentation.py          preprocessing, thresholding, front detection, alignment
  analysis.py              test boundaries, frame filtering, convergence, CDF, RZ width
  pipeline.py              the three stages: setup_run / subtract_background / extract_structure
  plots.py                 every figure + save_results / load_results / overlays
  export.py                per-frame PNG overlay writers, uint8 display conversion
  chemistry.py             OPTIONAL Cantera/SDToolbox parameter computation
run_analysis.py          main entry point — one dataset, full chain
run_background_subtract.py  background subtraction only
run_compare.py           overlay finished runs across datasets
```

Each dataset's results land in `OUTPUT_ROOT/<folder name>_Analysis/`.

---

## The three scripts

These are thin, editable wrappers. You configure a run by editing the `CONFIG`
(and, for `run_analysis.py`, `PHYSICAL`) block at the top of the file, then run it.
Every run archives a dated copy of the script and a dated `.txt` log alongside its
results, so a result can always be traced back to the code and parameters that
produced it.

### 1. `run_analysis.py` — full analysis of one dataset

The main workflow. Edit three blocks at the top:

- **CONFIG** — paths and run options:
  - `DATA_FOLDER` — folder holding the image sequence (frame `001` = background).
  - `OUTPUT_ROOT` — where the `<name>_Analysis` folder is written.
  - `SAVE_SUBTRACTED_IMAGES` — also write the background-subtracted frames as PNGs.
  - `EXPORT_FRAME_OVERLAYS` — write per-frame Front / Aligned / Structure PNGs
    (slow; a debugging aid, off by default in spirit).
  - `EXPORT_MATLAB_RESULTS` — also write `FrontData.mat` for MATLAB downstream code.
- **PHYSICAL** — the experiment's physical parameters (a `PhysicalParams`). Every
  field is required and must be `> 0`; the run stops with a clear error if a
  placeholder is left unfilled. Fields: `mixture`, `pressure_kpa`,
  `induction_length_m`, `u_cj`, `ea_rt`, `px2mm`, `cell_width_mm`.
- **EXTRACT** — algorithm knobs (an `ExtractionParams`). Defaults are usually fine;
  see *Key parameters* below.

Then:

```bash
uv run python run_analysis.py
```

The `main()` runs the pipeline (`setup_run → subtract_background →
extract_structure`), saves results, then calls the plot functions. Comment out
plot calls you don't want. Outputs:

- `Stats/Heatmap.png`, `ProbabilityCurve.png`, `Convergence.png`,
  `OtsuStability.png`, `CoverageMap.png`, `CDF.png`, `ReactionZoneWidth.png`
- `results.npz` + `results_summary.json` (used by `run_compare.py`)
- optional `FrontData.mat`, `Subtracted/`, and per-frame overlay folders
- the dated script archive and `.txt` summary log

### 2. `run_background_subtract.py` — background subtraction only

A standalone tool when you just want the background-subtracted, zero-median
normalized frames and none of the structural analysis. Subtraction and
normalization happen in float with no clamping; clamping only occurs when the
8-bit PNGs are written. Edit `DATA_FOLDER` and `OUTPUT_ROOT`, then:

```bash
uv run python run_background_subtract.py
```

Output: `OUTPUT_ROOT/<name>_Analysis/Subtracted/` plus provenance files.

### 3. `run_compare.py` — compare finished runs

Overlays curves and distributions from datasets you've **already analyzed**. It
loads each run from its `results.npz` + `results_summary.json` (via
`plots.load_results`) — no re-extraction and no original images needed. All
datasets are converted to a single common `SCHEME` so runs analyzed under
different schemes still compare cleanly.

Edit the `CONFIG` block:

- `RESULTS_FOLDERS` — list of `<name>_Analysis` folders to compare.
- `OUTPUT_FOLDER` — where the overlay figures go.
- `SCHEME` — common axis (`"cell_widths"` or `"induction_lengths"`).
- `CDF_XLIM`, `CDF_TARGET`, `PROB_DISTANCE_RANGE`, `RZ_CENSOR_MARGIN_PX`,
  `RZ_BINS`, `LABELS` — overlay-specific options.

```bash
uv run python run_compare.py
```

Produces `CDFOverlay.png`, `ProbabilityOverlay.png`, and
`ReactionZoneWidthOverlay.png` (each with a vector PDF sibling). Per-dataset
diagnostics (heatmap, coverage map, convergence) are *not* overlaid — they live
in each run's own `Stats/` folder.

---

## Input data conventions

- A dataset is a single folder of frames. Extensions are searched in priority
  order (`.bmp`, `.png`, `.tif`, `.tiff`); the first extension that matches wins,
  so mixed-format folders are never silently merged.
- Files are sorted by name. The trailing number in each filename is parsed as the
  frame number.
- **The frame whose trailing number is `1` ("frame 001") is the background.** If
  no such frame exists, the first file is used as background with a warning.
- Frames should be single-channel grayscale. An RGB fallback exists (it
  reproduces MATLAB `rgb2gray` Rec. 601 weights) but warns, since it shouldn't be
  needed.
- The folder name becomes the analysis name and the output subfolder name.

## Key parameters (`ExtractionParams`)

The defaults mirror the MATLAB version. The ones you're most likely to touch:

- `normalization_scheme` — `"cell_widths"` (default, x/λ) or
  `"induction_lengths"` (x/Δᵢ). **Every** distance in the run — the axis, the
  probe distances, plot ranges — is in these units.
- `probe_distances` — distances (in scheme units) at which probability is tracked
  for convergence, e.g. `(0.02, 0.05, 0.10)`.
- `threshold_mode`:
  - `"pooled"` (default) — one folded-Otsu threshold from the combined histogram
    of all frames. Best for homogeneous datasets.
  - `"per_frame"` — threshold recomputed per frame.
  - `"noise"` — threshold = `noise_k × σ_noise` (robust MAD estimate). This is
    contrast-*independent*, so "structure" means the same physical thing across
    datasets. **Recommended for cross-mixture comparison and faint (e.g.
    Ar-diluted) data** where Otsu's bimodality assumption breaks down. Use the
    *same* `noise_k` across datasets you intend to compare.
- `segmentation_source` — `"normalized"` (default; threshold the raw
  background-subtracted data) or `"enhanced"` (threshold after median filter +
  unsharp mask). Either way the threshold is computed and applied in the same
  space — the V5 bug (threshold from enhanced, applied to normalized) is fixed.
- `front_min_object_px` — drop connected components smaller than this from the
  *front-detection* mask only (suppresses faint specks that hijack the
  rightmost-pixel detector). Does not affect the structure statistics.
- `structure_min_object_px` — drop small components from the structure mask
  itself. **This changes the measured statistics by definition** — it's a
  deliberate redefinition of "structure" as features of ≥ N pixels. Important for
  the reaction-zone width metric in faint data. Use the same value across compared
  datasets.
- `front_range_fraction` — keep only frames whose median front lies in this FOV
  fraction (default `(0.20, 0.95)`).
- `min_coverage` — pixels with fewer real-data frames than this are masked from
  the heatmap.
- `rng_seed` — seed for the randomized accumulation order (reproducible
  convergence).

## Optional: computing physical parameters from chemistry

If you'd rather not type measured `u_cj` / `induction_length_m`, the optional
`chemistry.py` layer can compute them from a Cantera mechanism (requires Cantera
+ SDToolbox installed separately):

```python
from detrz_extraction import chemistry

u_cj = chemistry.compute_cj_speed(pressure_pa, 295.0, "H2:2 O2:1 N2:2", "gri30.yaml")
ind  = chemistry.compute_znd_induction_length(u_cj, pressure_pa, 295.0,
                                              "H2:2 O2:1 N2:2", "gri30.yaml")
```

Note `compute_cj_speed` takes pressure in **Pa**, not kPa.

## Typical first run

```bash
uv sync
# edit run_analysis.py: set DATA_FOLDER, OUTPUT_ROOT, and the PHYSICAL block
uv run python run_analysis.py
# inspect OUTPUT_ROOT/<name>_Analysis/Stats/*.png
# repeat for other datasets, then:
# edit run_compare.py: list the _Analysis folders
uv run python run_compare.py
```

## ALSO INCLUDED: reactionZone_v7.m

The matlab script used in the paper "Hydrodynamic Structure Statistics from Schlieren Imaging in Detonations of Varying Regularity"

## ALSO INCLUDED: example_dataset

An example dataset of schlieren images of a CH4_2O2 detonation at 11 kPa. 

