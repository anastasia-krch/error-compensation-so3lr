#!/usr/bin/env python3
"""Assemble the gas-phase dissociation scan from DFT, APPLE&P and AMS/SO3LR.

Consumes AMS/ADF, AMS/APPLE&P and AMS/SO3LR output produced by run_ams_singlepoint.sh
over the scan geometries. Every method is read from its AMS output at the same coordinates.

The conventional dE columns are referenced to the separated fragments computed with the
same method.  For SO3LR, the output also contains a pair-plateau-referenced curve.  Its
zero is the mean energy of the intact ion-pair geometries at 16 and 18 Angstrom, where the
available scan is flat.  This avoids comparing a distant, charge-transferred neutral pair
with separately evaluated fragments constrained to formal charges +1 and -1.

Manifest is the TSV run_ams_singlepoint.sh consumed, one row per geometry:

    xyz<TAB>charge<TAB>tag

Rows whose tag begins with the scan prefix are scan points, and the tag's numeric suffix is
the separation in Angstrom. The two fragment rows are matched by name.

Output is an npz consumed by import_charge_plot_data.py.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from charge_common import index_ams_logs

HARTREE_TO_KJ_MOL = 2625.499639479
EV_TO_KJ_MOL = 96.48533212
CATION_ATOMS = 19
PAIR_ATOMS = 34          # 19 (EMIM+) + 15 (NTf2-); guards against off-atom Lp sites
SCAN_PREFIX = "scan_d"
FRAGMENT_TAGS = ("frag_cation", "frag_anion")
PLATEAU_MIN_SEPARATION_A = 16.0


def apple_cation_charge(result) -> float:
    """Net cation charge from APPLE&P's own reported charges.

    APPLE&P monopoles are fixed by parameterisation, so this is constant across the scan;
    the length guard catches the off-atom lone-pair site padding the array, which would
    otherwise shift the slice.
    """
    if result is None or result.get("mulliken") is None:
        return np.nan
    charges = np.asarray(result["mulliken"])
    if len(charges) != PAIR_ATOMS:
        return np.nan
    return float(charges[:CATION_ATOMS].sum())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True,
                       help="TSV of xyz<TAB>charge<TAB>tag, as given to the AMS driver")
    parser.add_argument("--logdir", type=Path, required=True,
                       help="directory of AMS stdout logs carrying tag=/xc=/basis= lines")
    parser.add_argument("--dft-level", default="PBE0/TZ2P")
    parser.add_argument("--apple-level", default="APPLE&P/ff")
    parser.add_argument("--so3lr-level", default="SO3LR/ASE")
    parser.add_argument("--output", type=Path, required=True, help="npz to write")
    args = parser.parse_args()

    from ase.io import read

    index = index_ams_logs(args.logdir)
    if not index:
        raise SystemExit(f"no AMS logs with provenance lines found in {args.logdir}")
    print(f"indexed {len(index)} AMS results across "
          f"{len({level for _, level in index})} levels of theory")

    rows, fragments = [], {}
    for line in args.manifest.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        xyz, charge, tag = line.split("\t")
        atoms = read(args.manifest.parent / xyz if not Path(xyz).is_absolute() else xyz)
        so3lr = index.get((tag, args.so3lr_level))
        if so3lr is None:
            print(f"  {tag}: no {args.so3lr_level} result, skipped")
            continue

        if tag in FRAGMENT_TAGS:
            fragments[tag] = dict(
                so3lr=so3lr["energy_hartree"],
                dft=index.get((tag, args.dft_level), {}).get("energy_hartree", np.nan),
                apple=index.get((tag, args.apple_level), {}).get("energy_hartree", np.nan))
            continue
        if not tag.startswith(SCAN_PREFIX):
            continue

        dft = index.get((tag, args.dft_level))
        if dft is None:
            print(f"  {tag}: no {args.dft_level} result, skipped")
            continue
        apple = index.get((tag, args.apple_level))
        charges = so3lr["atomic_charges"]
        if charges is None or len(charges) != len(atoms):
            raise SystemExit(f"{tag}: SO3LR AMS charges do not match the {len(atoms)}-atom geometry")
        rows.append(dict(
            separation=float(tag[len(SCAN_PREFIX):]),
            e_dft=dft["energy_hartree"],
            e_apple=apple["energy_hartree"] if apple else np.nan,
            e_so3lr=so3lr["energy_hartree"],
            q_so3lr=float(charges[:CATION_ATOMS].sum()),
            q_mdcq=float(dft["mdc_q"][:CATION_ATOMS].sum())
            if dft["mdc_q"] is not None else np.nan,
            q_apple=apple_cation_charge(apple)))

    if not rows:
        raise SystemExit("no scan points with usable DFT results")
    missing = [tag for tag in FRAGMENT_TAGS if tag not in fragments]
    if missing:
        raise SystemExit(f"manifest has no fragment rows for {missing}; without them the "
                         "dissociation energies have no per-method zero")

    rows.sort(key=lambda row: row["separation"])
    reference = {key: sum(fragments[tag][key] for tag in FRAGMENT_TAGS)
                 for key in ("so3lr", "dft", "apple")}

    def relative(column, zero, scale):
        return (np.array([row[column] for row in rows]) - zero) * scale

    separation = np.array([row["separation"] for row in rows])
    e_so3lr = np.array([row["e_so3lr"] for row in rows])
    plateau_mask = separation >= PLATEAU_MIN_SEPARATION_A
    if plateau_mask.sum() < 2:
        raise SystemExit(
            f"need at least two scan points at or beyond {PLATEAU_MIN_SEPARATION_A:.1f} A "
            "to define the SO3LR pair-energy plateau")
    so3lr_plateau_hartree = float(e_so3lr[plateau_mask].mean())

    result = dict(
        r=separation,
        dE_pbe0=relative("e_dft", reference["dft"], HARTREE_TO_KJ_MOL),
        dE_ap=relative("e_apple", reference["apple"], HARTREE_TO_KJ_MOL),
        dE_ml=relative("e_so3lr", reference["so3lr"], HARTREE_TO_KJ_MOL),
        dE_ml_pair_plateau=(e_so3lr - so3lr_plateau_hartree) * HARTREE_TO_KJ_MOL,
        so3lr_pair_plateau_min_A=np.array(PLATEAU_MIN_SEPARATION_A),
        q_ml=np.array([row["q_so3lr"] for row in rows]),
        q_mdcq_pbe0=np.array([row["q_mdcq"] for row in rows]),
        q_ap=np.array([row["q_apple"] for row in rows]))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output, **result)
    print(f"wrote {args.output} ({len(rows)} scan points, "
          f"{result['r'].min():.1f}-{result['r'].max():.1f} A)")
    plateau_spread = np.ptp(result["dE_ml_pair_plateau"][plateau_mask])
    print(f"SO3LR pair-plateau zero: mean of {plateau_mask.sum()} points at "
          f"r >= {PLATEAU_MIN_SEPARATION_A:.1f} A "
          f"(spread {plateau_spread:.4f} kJ/mol)")

    # A fixed-charge force field cannot have a charged-fragment inconsistency, so its
    # curve must reach the Coulomb expectation -- which at the largest separation is NOT
    # zero. Printing both is what distinguishes a correct ionic dissociation limit from
    # SO3LR's flat plateau.
    if np.isfinite(result["dE_ap"]).any():
        ke = 14.399645351950548
        expected = -ke / result["r"].max() * EV_TO_KJ_MOL
        print(f"APPLE&P at {result['r'].max():.0f} A: {result['dE_ap'][-1]:+.1f} vs "
              f"-ke/r = {expected:+.1f} kJ/mol")
        spread = np.ptp(result["q_ap"][np.isfinite(result["q_ap"])])
        if spread > 1e-6:
            print(f"  WARNING: APPLE&P charge varies by {spread:.2e} e across the scan; "
                  "it is fixed by parameterisation, so check the atom typing")


if __name__ == "__main__":
    main()
