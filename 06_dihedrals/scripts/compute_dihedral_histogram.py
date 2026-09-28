#!/usr/bin/env python3
"""Compute a periodic NTf2 C-S-N-S / S-N-S-C histogram from reduced XYZ."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


N_ATOMS, N_ANIONS = 1250, 250


def dihedral(p0, p1, p2, p3, box):
    b0 = p1 - p0; b1 = p2 - p1; b2 = p3 - p2
    b0 -= box * np.rint(b0 / box); b1 -= box * np.rint(b1 / box); b2 -= box * np.rint(b2 / box)
    b1 /= np.linalg.norm(b1, axis=1)[:, None]
    v = b0 - (b0 * b1).sum(1)[:, None] * b1; w = b2 - (b2 * b1).sum(1)[:, None] * b1
    return np.degrees(np.arctan2((np.cross(b1, v) * w).sum(1), (v * w).sum(1))) % 360


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xyz", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bins", type=int, default=72)
    args = parser.parse_args()
    counts = np.zeros((args.bins, args.bins), dtype=np.int64); frames = 0
    with args.xyz.open() as source:
        while (line := source.readline()):
            if int(line) != N_ATOMS: raise ValueError("expected a 1,250-atom dihedral XYZ frame")
            comment = source.readline(); lattice = comment.split('Lattice="', 1)[1].split('"', 1)[0].split()
            box = np.array([float(lattice[0]), float(lattice[4]), float(lattice[8])])
            xyz = np.array([[float(value) for value in source.readline().split()[1:]] for _ in range(N_ATOMS)])
            first = dihedral(xyz[:250], xyz[250:500], xyz[500:750], xyz[750:1000], box)
            second = dihedral(xyz[250:500], xyz[500:750], xyz[750:1000], xyz[1000:1250], box)
            counts += np.histogram2d(first, second, bins=args.bins, range=((0, 360), (0, 360)))[0].astype(np.int64); frames += 1
    centres = (np.arange(args.bins) + .5) * 360 / args.bins
    table = pd.DataFrame([(x, y, int(counts[i, j]), counts[i, j] / counts.sum())
                          for i, x in enumerate(centres) for j, y in enumerate(centres)],
                         columns=["csns_deg", "snsc_deg", "count", "probability"])
    args.output.parent.mkdir(parents=True, exist_ok=True); table.to_csv(args.output, index=False)
    print(f"wrote {args.output} from {frames} frames")


if __name__ == "__main__": main()
