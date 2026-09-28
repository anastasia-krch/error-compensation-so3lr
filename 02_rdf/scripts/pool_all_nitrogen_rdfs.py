#!/usr/bin/env python3
"""Pool all 300 K LionAnalysis nitrogen RDFs into compact plot data.

Each RDF is pooled from raw counts, not from an unweighted average of
normalised seed curves. The source output tree is ``<root>/<ff>/300K/seed*/``
and force-field directories are apple, so3lr, and opls (CL&P).
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np

FORCE_FIELDS = {"APPLE&P": "apple", "SO3LR": "so3lr", "CL&P": "opls"}
PAIR_LABELS = {"anion_anion": "anion N--anion N", "cation_cation": "cation Nmid--cation Nmid",
               "cation_anion": "cation Nmid--anion N"}
STEPS = re.compile(r"nSteps\s*=\s*(\d+)")
DENSITY = re.compile(r"numberdens\.\s*=\s*([0-9.eE+-]+)")
N_IONS = 250


def read_curve(path: Path):
    header = path.read_text().splitlines()[0]
    steps, density = STEPS.search(header), DENSITY.search(header)
    if not steps or not density:
        raise ValueError(f"missing nSteps/density in {path}")
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 4], float(steps.group(1)), float(density.group(1))


def pool(files: list[Path]):
    if not files:
        raise FileNotFoundError("no source RDFs")
    radii = counts = None
    normalisation = 0.0
    for path in files:
        r, raw, steps, density = read_curve(path)
        if radii is None:
            radii, counts = r, np.zeros_like(raw)
        elif not np.array_equal(radii, r):
            raise ValueError(f"incompatible bins: {path}")
        counts += raw
        normalisation += steps * density
    dr = float(radii[1] - radii[0])
    shell = 4 * np.pi / 3 * ((radii + dr / 2) ** 3 - (radii - dr / 2) ** 3)
    return radii, counts / (N_IONS * shell * normalisation)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--temperature", type=int, default=300)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for force_field, directory in FORCE_FIELDS.items():
        for pair, label in PAIR_LABELS.items():
            files = sorted((args.root / directory / f"{args.temperature}K").glob(f"seed*/rdf_{pair}_*.txt"))
            r, g = pool(files)
            rows.extend({"force_field": force_field, "pair": pair, "distance_A": f"{x:.8f}", "g_r": f"{y:.8f}"}
                        for x, y in zip(r, g))
            print(f"{force_field:7s} {label:28s} {len(files)} seeds")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["force_field", "pair", "distance_A", "g_r"])
        writer.writeheader(); writer.writerows(rows)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
