#!/usr/bin/env python3
"""Export the four 310 K dihedral histograms used for the comparison figure."""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

N_BINS = 72
EDGES = np.linspace(0.0, 360.0, N_BINS + 1)
CENTRES = 0.5 * (EDGES[:-1] + EDGES[1:])


def read_histogram(path: Path) -> np.ndarray:
    histogram = np.load(path)
    if histogram.shape != (N_BINS, N_BINS) or histogram.sum() <= 0:
        raise ValueError(f"Expected a non-empty {N_BINS}x{N_BINS} histogram: {path}")
    return histogram.astype(float)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--apple-glob", required=True,
                        help="glob for one APPLE&P seed's histogram chunks")
    parser.add_argument("--so3lr", type=Path, required=True)
    parser.add_argument("--clp", type=Path, required=True)
    parser.add_argument("--clp-scaled", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    apple_files = [Path(path) for path in sorted(glob.glob(args.apple_glob))]
    if not apple_files:
        raise FileNotFoundError(args.apple_glob)
    histograms = {
        "APPLE&P": sum((read_histogram(path) for path in apple_files), np.zeros((N_BINS, N_BINS))),
        "SO3LR": read_histogram(args.so3lr),
        "CL&P": read_histogram(args.clp),
        "CL&P-scaled": read_histogram(args.clp_scaled),
    }
    phi_1, phi_2 = np.meshgrid(CENTRES, CENTRES, indexing="ij")
    rows = []
    bin_area = (EDGES[1] - EDGES[0]) ** 2
    for force_field, histogram in histograms.items():
        density = histogram / (histogram.sum() * bin_area)
        rows.append(pd.DataFrame({"force_field": force_field,
                                  "phi_1_deg": phi_1.ravel(),
                                  "phi_2_deg": phi_2.ravel(),
                                  "probability_density_deg_minus2": density.ravel()}))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(rows, ignore_index=True).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
