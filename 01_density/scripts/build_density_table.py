#!/usr/bin/env python3
"""Build the compact density plot table from archived NPT logs.

The input manifest is a local CSV with columns:
    force_field,engine,temperature_K,log,mass_amu

`log` paths are intentionally local and are never committed. The output has
only the values plotted in the paper: force field, temperature, density, SEM.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from extract_density import analyze


OUTPUT_FIELDS = [
    "force_field",
    "temperature_K",
    "density_g_cm3",
    "density_sem_g_cm3",
]
REQUIRED_FIELDS = {"force_field", "engine", "temperature_K", "log", "mass_amu"}


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = REQUIRED_FIELDS - fields
        if missing:
            raise ValueError(f"{path}: missing manifest columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path}: manifest has no simulation rows")
    return rows


def build(rows: list[dict[str, str]], fraction: float, blocks: int) -> list[dict[str, object]]:
    table = []
    for row in rows:
        log = Path(row["log"])
        if not log.is_file():
            raise FileNotFoundError(f"Archived log not found: {log}")
        engine = row["engine"].strip().lower()
        if engine not in {"ams", "so3lr"}:
            raise ValueError(f"{log}: engine must be 'ams' or 'so3lr', got {engine!r}")
        result = analyze(engine, log, float(row["mass_amu"]), fraction, blocks)
        table.append({
            "force_field": row["force_field"].strip(),
            "temperature_K": f"{float(row['temperature_K']):.0f}",
            "density_g_cm3": f"{result['density_mean_g_cm3']:.6f}",
            "density_sem_g_cm3": f"{result['density_sem_g_cm3']:.6f}",
        })
    return sorted(table, key=lambda item: (item["force_field"], float(item["temperature_K"])))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="local CSV listing archived NPT logs")
    parser.add_argument("--output", type=Path, required=True, help="plot-data CSV to create")
    parser.add_argument("--fraction", type=float, default=0.5, help="final trajectory fraction (default: 0.5)")
    parser.add_argument("--blocks", type=int, default=10, help="number of SEM blocks (default: 10)")
    args = parser.parse_args(argv)

    table = build(read_manifest(args.manifest), args.fraction, args.blocks)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(table)
    print(f"wrote {len(table)} plot-data rows to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
