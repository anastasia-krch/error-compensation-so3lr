#!/usr/bin/env python3
"""Assemble APPLE&P ion-pair XYZs from verified post-fix NVE chunks.


It creates a *working* source tree for LionAnalysis, never copies trajectories
into the publication repository.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path


def frames(path: Path):
    with path.open() as handle:
        while line := handle.readline():
            n_atoms = int(line)
            comment = handle.readline()
            body = [handle.readline() for _ in range(n_atoms)]
            if len(body) != n_atoms or not all(body):
                raise ValueError(f"truncated XYZ frame in {path}")
            if n_atoms != 500:
                raise ValueError(f"{path}: expected 500 atoms, found {n_atoms}")
            yield line, comment, body


def chunk_number(path: Path) -> int:
    match = re.search(r"chunk(\d+)_ionpair\.xyz$", path.name)
    if not match:
        raise ValueError(path)
    return int(match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help=".../nve_extracted directory")
    parser.add_argument("destination", type=Path, help="working LionAnalysis-root/apple directory")
    parser.add_argument("--temperature", type=int, default=300)
    parser.add_argument("--first-corrected-chunk", type=int, default=4)
    parser.add_argument("--seeds", type=int, nargs="+", default=(1, 2, 3, 4, 5))
    args = parser.parse_args()
    for seed in args.seeds:
        source_dir = args.source / f"{args.temperature}K_seed{seed}"
        chunk_paths = sorted(source_dir.glob("chunk*_ionpair.xyz"), key=chunk_number)
        chunk_paths = [path for path in chunk_paths if chunk_number(path) >= args.first_corrected_chunk]
        if not chunk_paths:
            raise FileNotFoundError(f"No chunks >= {args.first_corrected_chunk}: {source_dir}")
        expected = list(range(args.first_corrected_chunk, chunk_number(chunk_paths[-1]) + 1))
        found = [chunk_number(path) for path in chunk_paths]
        if found != expected:
            raise ValueError(f"Non-contiguous chunks for seed {seed}: {found}")
        output_dir = args.destination / f"{args.temperature}K" / f"seed{seed}"
        output_dir.mkdir(parents=True, exist_ok=True)
        output = output_dir / f"apple_{args.temperature}K_seed{seed}_ionpair.xyz"
        n_frames = 0
        with output.open("w") as handle:
            for i, path in enumerate(chunk_paths):
                for j, frame in enumerate(frames(path)):
                    if i and j == 0:
                        continue
                    handle.write(frame[0]); handle.write(frame[1]); handle.writelines(frame[2])
                    n_frames += 1
        print(f"seed {seed}: chunks {found[0]}–{found[-1]}, {n_frames} frames -> {output}")


if __name__ == "__main__":
    main()
