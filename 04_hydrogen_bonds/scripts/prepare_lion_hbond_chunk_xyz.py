#!/usr/bin/env python3
"""Stream selected one-based atom-ID ranges from a universal H-bond XYZ."""
from __future__ import annotations

import argparse
from pathlib import Path


UNIVERSAL_ATOMS = 4000


def atom_range(value: str) -> tuple[int, int]:
    first, last = (int(part) for part in value.split(":", 1))
    if not 1 <= first <= last <= UNIVERSAL_ATOMS:
        raise argparse.ArgumentTypeError("range must be within 1:4000")
    return first, last


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--range", type=atom_range, action="append", required=True,
                        help="inclusive atom-ID range, e.g. 751:1000; may be repeated")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    n_atoms = sum(last - first + 1 for first, last in args.range)
    frames = 0
    with args.input.open() as source, args.output.open("w") as destination:
        while (count_line := source.readline()):
            if int(count_line) != UNIVERSAL_ATOMS:
                raise ValueError(f"expected {UNIVERSAL_ATOMS} atoms, got {count_line.strip()}")
            comment = source.readline()
            coordinates = [source.readline() for _ in range(UNIVERSAL_ATOMS)]
            if any(not line for line in coordinates):
                raise ValueError("truncated XYZ frame")
            destination.write(f"{n_atoms}\n{comment}")
            for first, last in args.range:
                destination.writelines(coordinates[first - 1:last])
            frames += 1
    print(f"wrote {args.output} ({frames} frames, {n_atoms} atoms/frame)")


if __name__ == "__main__":
    main()
