#!/usr/bin/env python3
"""Run donor-resolved H-bond RDFs with LionAnalysis 0.5 safely.

LionAnalysis 0.5 silently reads just one frame from XYZ files above 500 atoms.
The universal reduced H-bond XYZ has 4,000 atoms, therefore this runner splits
each donor--acceptor calculation into 500-atom chunks, sums the raw pair counts,
and writes one merged Lion-style file per RDF. Chunk XYZ/config/log files are
temporary and removed after a successful merge.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


PAIRS = {
    "HCR_O": ((range(751, 1001),), (range(1501, 1751), range(1751, 2001), range(2001, 2251), range(2251, 2501)), 250, 1000),
    "HCR_F": ((range(751, 1001),), tuple(range(start, start + 250) for start in range(2501, 4001, 250)), 250, 1500),
    "HCW_O": ((range(1001, 1251), range(1251, 1501)), (range(1501, 1751), range(1751, 2001), range(2001, 2251), range(2251, 2501)), 500, 1000),
    "HCW_F": ((range(1001, 1251), range(1251, 1501)), tuple(range(start, start + 250) for start in range(2501, 4001, 250)), 500, 1500),
}


def specification(values: range) -> str:
    return f"{values.start}:{values.stop - 1}"


def parse_raw(path: Path) -> tuple[int, float, list[tuple[float, int]]]:
    header, *rows = path.read_text().splitlines()
    steps = int(re.search(r"nSteps = (\d+)", header).group(1))
    density = float(re.search(r"numberdens\. = ([0-9.eE+-]+)", header).group(1))
    values = []
    for row in rows:
        fields = row.split()
        values.append((float(fields[0]), int(round(float(fields[4])))))
    return steps, density, values


def merge(chunks: list[Path], output: Path, pair: str, n_donors: int, n_acceptors: int) -> None:
    parsed = [parse_raw(path) for path in chunks]
    steps = parsed[0][0]
    if any(item[0] != steps for item in parsed):
        raise ValueError("inconsistent frame counts between chunks")
    radii = [radius for radius, _ in parsed[0][2]]
    counts = [sum(item[2][index][1] for item in parsed) for index in range(len(radii))]
    density = sum(item[1] for item in parsed)
    donor, acceptor = pair.split("_")
    with output.open("w") as handle:
        handle.write("#r RDF iRDF PMF(kTunits) rawcount rawcountpertimestepperatom "
                     f"# from {donor} to {acceptor}; nSteps = {steps} size({donor}) = {n_donors} "
                     f"size({acceptor}) = {n_acceptors} numberdens. = {density:.12g}\n")
        for radius, count in zip(radii, counts):
            handle.write(f"{radius:.8g} 0 0 inf {count} {count / (steps * n_donors):.12g}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path, help="root containing <force-field>/300K/seed*/..._hbond.xyz")
    parser.add_argument("output_root", type=Path, help="external directory for merged Lion output")
    parser.add_argument("--temperature", type=int, default=300)
    parser.add_argument("--lion", type=Path, required=True)
    parser.add_argument("--force-field", action="append", choices=("so3lr", "clp", "apple", "clp-scaled"))
    parser.add_argument("--source-name", action="append", default=[], metavar="PUBLIC=ON_DISK",
                        help="e.g. clp=opls for input folder aliases")
    args = parser.parse_args()
    source_names = dict(item.split("=", 1) for item in args.source_name)
    selected = args.force_field or ("so3lr", "clp", "apple", "clp-scaled")
    script_dir = Path(__file__).parent
    prepare = script_dir / "prepare_lion_hbond_chunk_xyz.py"
    writer = script_dir / "write_lionanalysis_hbond_rdfs.py"
    for force_field in selected:
        source_dir = source_names.get(force_field, force_field)
        inputs = sorted((args.input_root / source_dir / f"{args.temperature}K").glob("seed*/*_hbond.xyz"))
        if not inputs:
            raise SystemExit(f"no input files for {force_field} in {args.input_root / source_dir}")
        for xyz in inputs:
            label = xyz.stem.removesuffix("_hbond")
            run_dir = args.output_root / force_field / f"{args.temperature}K" / xyz.parent.name
            run_dir.mkdir(parents=True, exist_ok=True)
            for pair, (donor_blocks, acceptor_blocks, n_donors, n_acceptors) in PAIRS.items():
                raw_chunks = []
                for donor_index, donor_range in enumerate(donor_blocks):
                    for acceptor_index, acceptor_range in enumerate(acceptor_blocks):
                        suffix = f"{pair}_d{donor_index}_a{acceptor_index}"
                        compact = run_dir / f"{suffix}.xyz"
                        config = run_dir / f"{suffix}.config"
                        raw = run_dir / f"rdf_{pair}_{label}_{suffix}.txt"
                        log = run_dir / f"{suffix}.out"
                        subprocess.run([sys.executable, str(prepare), "--input", str(xyz),
                                        "--range", specification(donor_range), "--range", specification(acceptor_range),
                                        "--output", str(compact)], check=True)
                        subprocess.run([sys.executable, str(writer), "--xyz", str(compact), "--label", f"{label}_{suffix}",
                                        "--pair", pair, "--n-donor", str(len(donor_range)),
                                        "--n-acceptor", str(len(acceptor_range)), "--output", str(config)], check=True)
                        with log.open("w") as handle:
                            subprocess.run([str(args.lion), str(config)], cwd=run_dir, stdout=handle,
                                           stderr=subprocess.STDOUT, check=True)
                        raw_chunks.append(raw)
                        compact.unlink()
                        config.unlink()
                    
                merge(raw_chunks, run_dir / f"rdf_{pair}_{label}.txt", pair, n_donors, n_acceptors)
                for raw in raw_chunks:
                    raw.unlink()
                print(f"finished {force_field}: {label} {pair}", flush=True)


if __name__ == "__main__":
    main()
