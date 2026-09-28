#!/usr/bin/env python3
"""Convert neutral-control references into AMS manifests before and after ADF optimisation."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def write_manifest(rows, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w") as handle:
        handle.write("xyz\tcharge\ttag\n")
        for xyz, charge, tag in rows:
            if not xyz.is_file():
                raise SystemExit(f"missing geometry for {tag}: {xyz}")
            handle.write(f"{xyz}\t{charge}\t{tag}\n")
    print(f"wrote {destination} ({len(rows)} geometries)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True,
                        help="directory holding ref_manifest.tsv and ref_<species>.xyz")
    parser.add_argument("--initial-output", type=Path, required=True,
                        help="AMS manifest for ADF geometry optimisation")
    parser.add_argument("--optimised-root", type=Path,
                        help="ADF run root containing <species>/ams.results/output.xyz")
    parser.add_argument("--so3lr-output", type=Path,
                        help="AMS manifest for SO3LR single points at the optimised geometries")
    args = parser.parse_args()
    if (args.optimised_root is None) != (args.so3lr_output is None):
        raise SystemExit("--optimised-root and --so3lr-output must be provided together")

    reference = args.reference_dir / "ref_manifest.tsv"
    with reference.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if not records or set(records[0]) < {"name", "charge"}:
        raise SystemExit(f"{reference} must have name and charge columns")

    initial = [(args.reference_dir / f"ref_{row['name']}.xyz", row["charge"], row["name"])
               for row in records]
    write_manifest(initial, args.initial_output)
    if args.optimised_root:
        optimised = [(args.optimised_root / row["name"] / "ams.results" / "output.xyz",
                      row["charge"], row["name"])
                     for row in records]
        write_manifest(optimised, args.so3lr_output)


if __name__ == "__main__":
    main()
