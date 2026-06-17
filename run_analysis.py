#!/usr/bin/env python3
"""run_analysis.py — wavefront detection, alignment, and reaction zone analysis.

Author: Raj Punna / Claude
Date: June 2026
Port of reactionZone_v5.m.

FIRST RUN: fill in the CONFIG variables and PHYSICAL parameters below
(EXTRACT defaults are usually fine), then:  uv run python run_analysis.py

The script is three stages — comment out what you don't need:
  1. subtract_background  — standalone tool; set SAVE_SUBTRACTED_IMAGES to
                            keep the background-subtracted dataset
  2. extract_structure    — the full front-detection / statistics chain
  3. plots.*              — pick which analyses/figures you want

OUTPUT:
  [OUTPUT_ROOT]/[dataset name]_Analysis/
  ├── Subtracted/                    <- background-subtracted frames (opt-in)
  ├── Fronts/ Aligned/ Structure/    <- per-frame overlays (opt-in)
  ├── Stats/                         <- the figures you asked for
  ├── results.npz / results_summary.json / FrontData.mat (opt-in)
  ├── run_analysis_YYYY-MM-DD.py     <- script archive (provenance)
  └── run_analysis_YYYY-MM-DD.txt  <- summary log
"""

from pathlib import Path

from detrz_extraction.parameters import ExtractionParams, PhysicalParams
from detrz_extraction import pipeline, plots

# ========== CONFIG —    paths and run options (edit me) ==========

DATA_FOLDER = Path('/Users/19rwrp/Desktop/Project/Opensource/datasets/CH4_2O2 11')  # folder with the image sequence
OUTPUT_ROOT = Path("output_final")                 # results: output/<folder name>_Analysis
SAVE_SUBTRACTED_IMAGES = False               # keep the background-subtracted dataset
EXPORT_FRAME_OVERLAYS = True                # per-frame PNGs: slow; enable for debugging
EXPORT_MATLAB_RESULTS = False                # write FrontData.mat for MATLAB downstream

# ========== PHYSICAL — test parameters for THIS dataset (edit me) ==========
# Every field is required; the run stops with a clear error if any is left at 0.

PHYSICAL = PhysicalParams(
    mixture="CH4_2O2",       # mixture label, e.g. "2H2_O2_2N2"
    pressure_kpa=11,          # initial pressure [kPa]
    induction_length_m=0.003583,  # ZND induction length [m]
    u_cj=2291,                # CJ velocity [m/s]
    ea_rt=12.33,                  # reduced activation energy Ea/RT
    px2mm=0.2497,                  # pixel scale [mm/px]
    cell_width_mm=64.3,         # detonation cell width [mm]
)

# ========== EXTRACT — algorithm knobs (defaults usually fine) ==========

EXTRACT = ExtractionParams()
# All distances (probes, plot ranges, dist_range) share ONE normalization
# scheme: "cell_widths" (default, x/lambda) or "induction_lengths" (x/Delta_i).
# Override individual fields like so:
# EXTRACT.normalization_scheme = "induction_lengths"
# EXTRACT.probe_distances = (0.02, 0.05, 0.10, 0.20)   # in scheme units
# EXTRACT.front_range_fraction = (0.30, 0.70)          # old conservative window
# For faint-schlieren datasets / cross-mixture comparisons (see CHANGELOG #21):
#EXTRACT.threshold_mode = "noise"        # contrast-independent k*sigma criterion
#EXTRACT.noise_k = 4.0                   # SAME k across compared datasets
#EXTRACT.front_min_object_px = 8         # noise ahead of front hijacks detection
#EXTRACT.structure_min_object_px = 5     # redefines structure as >=N px features

# -- OPTIONAL chemistry layer (requires Cantera + SDToolbox; see chemistry.py) --
# Computes the PHYSICAL values above from a mechanism instead of typing them:
# from detrz_extraction import chemistry
# u_cj = chemistry.compute_cj_speed(pressure_pa, 295.0, "H2:2 O2:1 N2:2", "gri30.yaml")
# induction_length_m = chemistry.compute_znd_induction_length(u_cj, ...)


def main() -> None:
    run = pipeline.setup_run(DATA_FOLDER, OUTPUT_ROOT, __file__)

    # --- Stage 1: background subtraction (usable standalone) ---
    dataset = pipeline.subtract_background(
        DATA_FOLDER,
        extensions=EXTRACT.frame_extensions,
        save_images=SAVE_SUBTRACTED_IMAGES,
        run=run,
    )

    # --- Stage 2: structural extraction + statistics ---
    results = pipeline.extract_structure(
        dataset, EXTRACT, PHYSICAL, run,
        export_frame_overlays=EXPORT_FRAME_OVERLAYS,
    )
    plots.save_results(results, run, matlab=EXPORT_MATLAB_RESULTS)

    # --- Stage 3: pick your analyses ---
    plots.plot_heatmap(results, run)
    plots.plot_probability_curve(results, run)
    # plots.plot_probability_curve(results, run, distance_range=(0, 0.5),  # scheme units
    #                              filename="ProbabilityCurve_near_front.png")
    plots.plot_convergence(results, run)
    plots.plot_otsu_stability(results, run)
    plots.plot_coverage_map(results, run)
    plots.plot_cdf_overlay([results], run.folders["stats"] / "CDF.png")
    plots.plot_reaction_zone_width(results, run)
    # Visualize exactly which pixels feed that distribution (per-frame PNGs):
    plots.export_reaction_zone_pixels(dataset, results, run, montage=True)

    # Structural pixels colored by distance from the front (slow; per-frame PNGs):
    plots.export_distance_colored_frames(dataset, results, run,
                                          dist_range=(0, 1.2), montage=True)

    # To compare against previously analyzed datasets, see run_compare.py

    plots.log_summary(results, run)


if __name__ == "__main__":
    main()
