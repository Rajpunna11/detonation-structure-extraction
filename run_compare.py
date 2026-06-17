#!/usr/bin/env python3
"""run_compare.py — compare results across previously analyzed datasets.

Author: Raj Punna / Claude
Date: June 2026

Loads finished runs from their analysis folders (results.npz +
results_summary.json — no re-extraction, no original images needed) and
overlays them on common axes. All datasets are converted to SCHEME units
regardless of the scheme each was analyzed with, so mixed runs compare
cleanly.

PRINCIPLE: every analysis that outputs a curve or distribution has an
*_overlay function in detrz_extraction.plots taking a list of results
objects (in-memory ExtractionResults and/or load_results() bundles):
    plot_cdf_overlay                  - CDF of P(x)
    plot_probability_overlay          - P(x) vs distance
    plot_reaction_zone_width_overlay  - reaction-zone width distributions
Per-dataset images (heatmap, coverage map, convergence) are diagnostics of
a single run and live in that run's Stats/ folder instead.

FIRST RUN: list the analysis folders to compare, then:
    uv run python run_compare.py
"""

from pathlib import Path

from detrz_extraction import plots

# ========== CONFIG (edit me) ==========

RESULTS_FOLDERS = [
    Path('/Users/19rwrp/Desktop/Project/Opensource/detrz-extraction/output_noise/CH4_2O2 11_Analysis'),
    Path('/Users/19rwrp/Desktop/Project/Opensource/detrz-extraction/output_noise/2H2_O2_3Ar 6_Analysis'),
    Path('/Users/19rwrp/Desktop/Project/Opensource/detrz-extraction/output_noise/2H2_O2_2N2 10_Analysis'),
]
OUTPUT_FOLDER = Path("output_noise/Comparison")

SCHEME = "cell_widths"        # common axis: "cell_widths" or "induction_lengths"
CDF_XLIM = (0.0, 0.75)        # CDF x-axis window (scheme units)
CDF_TARGET = 0.98             # CDF reference line
PROB_DISTANCE_RANGE = (0,0.2)    # e.g. (0, 0.5) to zoom P(x); None = full axis
RZ_CENSOR_MARGIN_PX = 2       # censoring margin for the width distributions
RZ_BINS = 60                  # histogram bins for the width overlay
LABELS = None                 # e.g. ["20 kPa", "40 kPa"]; None = dataset names


def main() -> None:
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

    bundles = [plots.load_results(folder) for folder in RESULTS_FOLDERS]
    print(f"Loaded {len(bundles)} datasets: {[b.dataset_name for b in bundles]}")

    p1 = plots.plot_cdf_overlay(bundles, OUTPUT_FOLDER / "CDFOverlay.png",
                                scheme=SCHEME, cdf_target=CDF_TARGET,
                                xlim=CDF_XLIM, labels=LABELS)
    p2 = plots.plot_probability_overlay(bundles, OUTPUT_FOLDER / "ProbabilityOverlay.png",
                                        scheme=SCHEME,
                                        distance_range=PROB_DISTANCE_RANGE, labels=LABELS)
    p3, rz_stats = plots.plot_reaction_zone_width_overlay(
        bundles, OUTPUT_FOLDER / "ReactionZoneWidthOverlay.png",
        scheme=SCHEME, censor_margin_px=RZ_CENSOR_MARGIN_PX, bins=RZ_BINS, labels=LABELS)
    for b, s in zip(bundles, rz_stats):
        print(f"  {b.dataset_name}: N={s['widths'].size}, "
              f"censored {s['censored_fraction'] * 100:.1f}%")
    print(f"Saved: {p1}\nSaved: {p2}\nSaved: {p3}")
    print("(Vector PDF siblings saved alongside each PNG.)")


if __name__ == "__main__":
    main()
