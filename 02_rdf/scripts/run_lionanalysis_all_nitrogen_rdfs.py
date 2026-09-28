#!/usr/bin/env python3
"""Run all three nitrogen RDFs for every reduced 300 K trajectory.

The runner writes its config and LionAnalysis outputs next to each source
reduced XYZ, outside the publication repository. Use
``pool_all_nitrogen_rdfs.py`` afterwards to create compact plot data.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


FORCE_FIELDS = ("apple", "so3lr", "opls")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="Lion reduced-trajectory root")
    parser.add_argument("--temperature", type=int, default=300)
    parser.add_argument("--lion", type=Path, required=True, help="lionanalysis.exe")
    parser.add_argument("--writer", type=Path, default=Path(__file__).with_name("write_lionanalysis_all_nitrogen_rdfs.py"))
    parser.add_argument("--force-field", choices=FORCE_FIELDS, action="append")
    args = parser.parse_args()
    force_fields = args.force_field or FORCE_FIELDS
    jobs = []
    for ff in force_fields:
        for xyz in sorted((args.root / ff / f"{args.temperature}K").glob("seed*/*_ionpair.xyz")):
            label = xyz.stem.removesuffix("_ionpair")
            jobs.append((ff, xyz, label))
    if not jobs:
        raise SystemExit("no reduced ion-pair XYZ files found")
    for number, (ff, xyz, label) in enumerate(jobs, 1):
        directory = xyz.parent
        config = directory / f"all_nitrogen_rdfs_{label}.config"
        print(f"[{number}/{len(jobs)}] {ff}: {label}", flush=True)
        subprocess.run([sys.executable, str(args.writer), "--xyz", str(xyz), "--label", label, "--output", str(config)], check=True)
        with (directory / f"all_nitrogen_rdfs_{label}.out").open("w") as log:
            subprocess.run([str(args.lion), str(config)], cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True)


if __name__ == "__main__":
    main()
