#!/usr/bin/env python3
"""Molecular dipoles for the neutral-analogue control: AMS/SO3LR against DFT.

Consumes AMS/ADF and AMS/SO3LR output produced by run_ams_singlepoint.sh (geometry
optimisations for the neutral analogues, single points for the ion-pair snapshots). Both
methods are evaluated at the same coordinates and are read from AMS logs.

Two npz files are written, both consumed by import_charge_plot_data.py:

    <out>/t1_neutral_control.npz   <species>_dipoles = [mu_DFT, mu_SO3LR]
    <out>/pair_comparison.npz      SO3LR_geoms, columns
                                   [q_MDCq, q_Mulliken, q_SO3LR, mu_DFT, mu_SO3LR]

The dipole is the comparison of choice because it needs no partitioning of the electron
density.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from charge_common import index_ams_logs

# The seven neutral analogues carry the same functional groups as the ions.
NEUTRALS = ("ch3cn", "ch3cooh", "methylimidazole", "ch3sh", "ch3sch3", "h2co", "hntf2")
CATION_ATOMS = 19       # the ion-pair xyz files list the cation first


def optimised_geometry(root: Path, species: str):
    """The structure AMS writes beside the run, for a species that was optimised."""
    from ase.io import read

    path = root / species / "ams.results" / "output.xyz"
    if not path.is_file():
        raise SystemExit(f"missing optimised geometry: {path}")
    return read(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry-dir", type=Path, required=True,
                       help="AMS run root holding <species>/ams.results/output.xyz")
    parser.add_argument("--logdir", type=Path, required=True,
                       help="directory of AMS stdout logs carrying tag=/xc=/basis= lines")
    parser.add_argument("--pair-geometry", type=Path, nargs="*", default=[],
                       help="ion-pair snapshot XYZ files, tagged by their file stem")
    parser.add_argument("--level", default="PBE/TZP",
                       help="level of theory to read from the log index")
    parser.add_argument("--so3lr-level", default="SO3LR/ASE")
    parser.add_argument("--out", type=Path, required=True, help="directory for the npz")
    args = parser.parse_args()

    index = index_ams_logs(args.logdir)
    if not index:
        raise SystemExit(f"no AMS logs with provenance lines found in {args.logdir}")
    print(f"indexed {len(index)} AMS results from {args.logdir}")

    args.out.mkdir(parents=True, exist_ok=True)

    neutral: dict[str, np.ndarray] = {}
    for species in NEUTRALS:
        result = index.get((species, args.level))
        so3lr = index.get((species, args.so3lr_level))
        if result is None or so3lr is None:
            print(f"  {species:18s} NO USABLE AMS RESULT at both levels -- skipped")
            continue
        optimised_geometry(args.geometry_dir, species)  # refuse a log/geometry mismatch
        mu_so3lr = so3lr["dipole_D"]
        neutral[f"{species}_dipoles"] = np.array([result["dipole_D"], mu_so3lr])
        if result["mdc_q"] is not None:
            neutral[f"{species}_dft_mdc_q"] = result["mdc_q"]
        if so3lr["atomic_charges"] is None:
            raise SystemExit(f"{species}: SO3LR AMS result has no atomic charges")
        neutral[f"{species}_so3lr_q"] = so3lr["atomic_charges"]
        print(f"  {species:18s} mu_DFT {result['dipole_D']:8.3f}  "
              f"mu_SO3LR {mu_so3lr:8.3f}  ratio {mu_so3lr / result['dipole_D']:.3f}")

    path = args.out / "t1_neutral_control.npz"
    np.savez(path, **neutral)
    print(f"wrote {path}")

    if args.pair_geometry:
        from ase.io import read

        rows = []
        for geometry in args.pair_geometry:
            tag = geometry.stem
            result = index.get((tag, args.level))
            so3lr = index.get((tag, args.so3lr_level))
            if result is None or so3lr is None:
                print(f"  {tag:18s} NO USABLE AMS RESULT at both levels -- skipped")
                continue
            atoms = read(geometry)
            charges = so3lr["atomic_charges"]
            if charges is None or len(charges) != len(atoms):
                raise SystemExit(f"{tag}: SO3LR AMS charges do not match the {len(atoms)}-atom geometry")
            mu_so3lr = so3lr["dipole_D"]
            # Net cation charge under each scheme, then the two partition-free dipoles.
            rows.append([
                float(result["mdc_q"][:CATION_ATOMS].sum())
                if result["mdc_q"] is not None else np.nan,
                float(result["mulliken"][:CATION_ATOMS].sum())
                if result["mulliken"] is not None else np.nan,
                float(charges[:CATION_ATOMS].sum()),
                float(result["dipole_D"]), mu_so3lr])
            print(f"  {tag:18s} q_SO3LR {rows[-1][2]:+.3f}  "
                  f"mu_DFT {rows[-1][3]:7.3f}  mu_SO3LR {mu_so3lr:7.3f}")
        path = args.out / "pair_comparison.npz"
        np.savez(path, SO3LR_geoms=np.array(rows))
        print(f"wrote {path} ({len(rows)} pairs)")


if __name__ == "__main__":
    main()
