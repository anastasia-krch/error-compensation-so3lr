#!/usr/bin/env python3
"""Write one LionAnalysis input for all N(midpoint)/anion-N RDFs.

The input ``*_nmid_anionN.xyz`` has anion central N atoms 1--250 and cation
ring-N midpoint atoms 251--500. Three RDFs are written from exactly those
coordinates: anion--anion, cation--cation, and cation--anion.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def ids(first: int, last: int) -> str:
    return " ".join(map(str, range(first, last + 1)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xyz", type=Path, required=True)
    parser.add_argument("--label", required=True, help="e.g. so3lr_300K_seed0")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-distance-A", type=float, default=15.0)
    parser.add_argument("--resolution-A", type=float, default=0.02)
    args = parser.parse_args()
    prefix = f"rdf_{{pair}}_{args.label}.txt"
    config = "\n".join([
        f"DumpFile {args.xyz}", "DumpFileFormat asexyz", "Overwrite", "Threads 1", "CoutFrequency 2000",
        f"DefineGroup AnionN LIST {ids(1, 250)}",
        f"DefineGroup CationNmid LIST {ids(251, 500)}",
        f"RDF AnionN AnionN {prefix.format(pair='anion_anion')} MaxDist {args.max_distance_A:.4f} Resolution {args.resolution_A:.4f}",
        f"RDF CationNmid CationNmid {prefix.format(pair='cation_cation')} MaxDist {args.max_distance_A:.4f} Resolution {args.resolution_A:.4f}",
        f"RDF AnionN CationNmid {prefix.format(pair='cation_anion')} MaxDist {args.max_distance_A:.4f} Resolution {args.resolution_A:.4f}",
        "",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(config)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
