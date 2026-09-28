#!/usr/bin/env python3
"""Import compact 300 K CL&P-scaled nitrogen RDFs into the paper plot table.

The source directory is an external analysis result. Only the three plotted
curves are copied into the repository; no trajectory or binary intermediate is
imported. The curves must use the documented direct reference-centred
normalisation and 0.02 A bin centres.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

FILES = {
    "anion_anion": "rdf_anion_anion_300K_seed0.txt",
    "cation_cation": "rdf_cation_cation_300K_seed0.txt",
    "cation_anion": "rdf_cation_anion_300K_seed0.txt",
}
LABEL = "CL&P-scaled"


def read_curve(path: Path) -> np.ndarray:
    values = np.loadtxt(path, comments="#", ndmin=2)
    if values.shape != (751, 2):
        raise ValueError(f"{path}: expected 751 distance/g(r) rows, found {values.shape}")
    expected = 0.01 + 0.02 * np.arange(751)
    if not np.allclose(values[:, 0], expected, rtol=0.0, atol=1e-8):
        raise ValueError(f"{path}: expected 0.02 A bin centres from 0.01 to 15.01 A")
    if not np.all(np.isfinite(values)) or np.any(values[:, 1] < 0.0):
        raise ValueError(f"{path}: RDF must be finite and non-negative")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path, help="external 300K_seed0 scaled-analysis directory")
    parser.add_argument("--output", type=Path, default=Path("../data/nitrogen_rdfs_300K.csv"))
    args = parser.parse_args()
    curves = {pair: read_curve(args.source_dir / name) for pair, name in FILES.items()}
    existing: list[dict[str, str]] = []
    if args.output.exists():
        with args.output.open(newline="") as handle:
            existing = [row for row in csv.DictReader(handle) if row["force_field"] != LABEL]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("force_field", "pair", "distance_A", "g_r"))
        writer.writeheader()
        writer.writerows(existing)
        for pair, values in curves.items():
            for distance, rdf in values:
                writer.writerow({"force_field": LABEL, "pair": pair,
                                 "distance_A": f"{distance:.8f}", "g_r": f"{rdf:.8f}"})
    print(f"Imported {len(FILES)} CL&P-scaled curves into {args.output}")


if __name__ == "__main__":
    main()
