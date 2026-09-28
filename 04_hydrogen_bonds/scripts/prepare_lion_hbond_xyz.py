#!/usr/bin/env python3
"""Stream one compact H-bond donor--acceptor XYZ from the universal reduced XYZ.

The universal file contains 4,000 atoms per frame. LionAnalysis 0.5 does not
reliably iterate this size, so this script produces one temporary input for one
RDF pair: HCR--O, HCR--F, HCW--O, or HCW--F. It preserves every sampled frame
and its lattice line. The temporary file is safe to delete after Lion finishes.
"""
from __future__ import annotations

import argparse
from pathlib import Path


UNIVERSAL_ATOMS = 4000
PAIR_SLICES = {
    "HCR_O": ((750, 1000), (1500, 2500)),
    "HCR_F": ((750, 1000), (2500, 4000)),
    "HCW_O": ((1000, 1500), (1500, 2500)),
    "HCW_F": ((1000, 1500), (2500, 4000)),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--pair", choices=PAIR_SLICES, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    first, second = PAIR_SLICES[args.pair]
    n_atoms = (first[1] - first[0]) + (second[1] - second[0])
    frames = 0
    with args.input.open() as source, args.output.open("w") as destination:
        while True:
            count_line = source.readline()
            if not count_line:
                break
            if int(count_line) != UNIVERSAL_ATOMS:
                raise ValueError(f"expected {UNIVERSAL_ATOMS} atoms, got {count_line.strip()}")
            comment = source.readline()
            coordinates = [source.readline() for _ in range(UNIVERSAL_ATOMS)]
            if any(not line for line in coordinates):
                raise ValueError("truncated XYZ frame")
            destination.write(f"{n_atoms}\n{comment}")
            destination.writelines(coordinates[first[0]:first[1]])
            destination.writelines(coordinates[second[0]:second[1]])
            frames += 1
    print(f"wrote {args.output} ({frames} frames, {n_atoms} atoms/frame)")


if __name__ == "__main__":
    main()
