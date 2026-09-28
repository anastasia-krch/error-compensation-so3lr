#!/usr/bin/env python3
"""Pool donor-resolved LionAnalysis H-bond RDF raw counts into plot data."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np


FORCE_FIELDS = {"APPLE&P": "apple", "SO3LR": "so3lr", "CL&P": "clp", "CL&P-scaled": "clp-scaled"}
PAIRS = {
    "HCR_O": (250, 1000), "HCR_F": (250, 1500),
    "HCW_O": (500, 1000), "HCW_F": (500, 1500),
}
STEPS = re.compile(r"nSteps\s*=\s*(\d+)")
DENSITY = re.compile(r"numberdens\.\s*=\s*([0-9.eE+-]+)")


def read_curve(path: Path) -> tuple[np.ndarray, np.ndarray, float, float]:
    header = path.read_text().splitlines()[0]
    steps, density = STEPS.search(header), DENSITY.search(header)
    if not steps or not density:
        raise ValueError(f"missing nSteps/numberdens in {path}")
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 4], float(steps.group(1)), float(density.group(1))


def pool(files: list[Path], n_donors: int) -> tuple[np.ndarray, np.ndarray]:
    if not files:
        raise FileNotFoundError("no source RDFs")
    radii = counts = None
    normalisation = 0.0
    for path in files:
        radius, raw, steps, density = read_curve(path)
        if radii is None:
            radii, counts = radius, np.zeros_like(raw)
        elif not np.array_equal(radii, radius):
            raise ValueError(f"incompatible bins: {path}")
        counts += raw
        normalisation += steps * density
    dr = float(radii[1] - radii[0])
    shell = 4 * np.pi / 3 * ((radii + dr / 2) ** 3 - (radii - dr / 2) ** 3)
    return radii, counts / (n_donors * shell * normalisation)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--temperature", type=int, default=300)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for force_field, directory in FORCE_FIELDS.items():
        for pair, (n_donors, _) in PAIRS.items():
            files = sorted((args.root / directory / f"{args.temperature}K").glob(f"seed*/rdf_{pair}_*.txt"))
            radii, rdf = pool(files, n_donors)
            rows.extend({"force_field": force_field, "pair": pair,
                         "distance_A": f"{r:.8f}", "g_r": f"{g:.8f}"}
                        for r, g in zip(radii, rdf))
            print(f"{force_field:12s} {pair:5s} {len(files)} seeds")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["force_field", "pair", "distance_A", "g_r"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
