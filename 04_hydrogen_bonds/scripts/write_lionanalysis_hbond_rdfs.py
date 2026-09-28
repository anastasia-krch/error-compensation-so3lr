#!/usr/bin/env python3
"""Write one LionAnalysis configuration for a compact donor--acceptor RDF XYZ."""
from __future__ import annotations

import argparse
from pathlib import Path


PAIR_GROUPS = {
    "HCR_O": ("HCR", 250, "Oacc", 1000),
    "HCR_F": ("HCR", 250, "Facc", 1500),
    "HCW_O": ("HCW", 500, "Oacc", 1000),
    "HCW_F": ("HCW", 500, "Facc", 1500),
}


def ids(first: int, last: int) -> str:
    return " ".join(map(str, range(first, last + 1)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xyz", type=Path, required=True)
    parser.add_argument("--label", required=True, help="e.g. clp_300K_seed0")
    parser.add_argument("--pair", choices=PAIR_GROUPS, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n-donor", type=int, help="donors present in compact XYZ")
    parser.add_argument("--n-acceptor", type=int, help="acceptors present in compact XYZ")
    parser.add_argument("--max-distance-A", type=float, default=8.0)
    parser.add_argument("--resolution-A", type=float, default=0.02)
    args = parser.parse_args()
    donor_name, default_donors, acceptor_name, default_acceptors = PAIR_GROUPS[args.pair]
    n_donor = args.n_donor or default_donors
    n_acceptor = args.n_acceptor or default_acceptors
    if n_donor < 1 or n_acceptor < 1:
        raise SystemExit("donor and acceptor counts must be positive")
    config = "\n".join([
        f"DumpFile {args.xyz}", "DumpFileFormat asexyz", "Overwrite", "Threads 1", "CoutFrequency 2000",
        f"DefineGroup {donor_name} LIST {ids(1, n_donor)}",
        f"DefineGroup {acceptor_name} LIST {ids(n_donor + 1, n_donor + n_acceptor)}",
        f"RDF {donor_name} {acceptor_name} rdf_{args.pair}_{args.label}.txt "
        f"MaxDist {args.max_distance_A:.4f} Resolution {args.resolution_A:.4f}",
        "",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(config)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
