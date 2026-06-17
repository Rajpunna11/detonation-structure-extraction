#!/usr/bin/env python3
"""run_background_subtract.py — standalone background subtraction.

Author: Raj Punna / Claude
Date: June 2026

Subtracts the background frame (frame 001 by convention) from every frame in
a sequence and zero-median normalizes — no structural extraction. Direct
subtraction + normalization in float, clamped only at PNG export.

FIRST RUN: set DATA_FOLDER and OUTPUT_ROOT, then:
    uv run python run_background_subtract.py

OUTPUT:
  [OUTPUT_ROOT]/[dataset name]_Analysis/
  ├── Subtracted/                          <- 8-bit PNGs of every frame
  ├── run_background_subtract_*.py / .txt  <- provenance + log
"""

from pathlib import Path

from detrz_extraction import pipeline

# ========== CONFIG (edit me) ==========

DATA_FOLDER = Path("examples/demo_dataset")  # folder with the image sequence
OUTPUT_ROOT = Path("output")                 # results: output/<folder name>_Analysis


def main() -> None:
    run = pipeline.setup_run(DATA_FOLDER, OUTPUT_ROOT, __file__)
    dataset = pipeline.subtract_background(DATA_FOLDER, save_images=True, run=run)
    run.log.info("Done: %d frames subtracted (%s background).",
                 dataset.n_frames, dataset.bg_name)


if __name__ == "__main__":
    main()
