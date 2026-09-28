#!/usr/bin/env python3
"""Flatten the charge-analysis caches into the committed CSV plot data."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# The npz caches written by the three compute_*.py scripts. They are binary and large
# relative to the numbers actually plotted, and *.npz is gitignored repo-wide, so the
# committed artefacts are the CSVs this script emits.
CONVERGENCE_NPZ = "charge_mechanism_300K.npz"
NEUTRAL_NPZ = "t1_neutral_control.npz"
PAIR_NPZ = "pair_comparison.npz"
PES_NPZ = "pes_scan_cache.npz"

# Seven neutral analogues, in the order the figure lists them, with the DFT-optimised
# geometry used for the structure thumbnail.
NEUTRALS = [
    ("ch3cn", "CH$_3$CN"),
    ("ch3cooh", "CH$_3$COOH"),
    ("methylimidazole", "1-methylimidazole"),
    ("ch3sh", "CH$_3$SH"),
    ("ch3sch3", "CH$_3$SCH$_3$"),
    ("h2co", "H$_2$CO"),
    ("hntf2", "HNTf$_2$"),
]
PAIR_SPECIES = "emim_ntf2_pair"
PAIR_LABEL = "EMIM$^+$/NTf$_2^-$"
PAIR_GEOMETRY = "pair_nve_1_c205_a208.xyz"

RECEPTIVE_FIELD_A = 13.5    # num_layers (3) x cutoff (4.5 A)


def _header(path: Path, lines: list[str]) -> None:
    """Prepend a `#` provenance preamble, read back with pd.read_csv(comment='#')."""
    body = path.read_text()
    path.write_text("".join(f"# {line}\n" for line in lines) + body)


def convergence_table(source: Path) -> pd.DataFrame:
    """Per-cluster-size mean and standard error of |q_ion|, plus the bulk reference.

    The mean/SEM reduction over the six central ion-pair environments is done here rather
    than in the notebook, so the committed CSV is what the figure plots and the notebook
    stays a plotter.
    """
    d = np.load(source / CONVERGENCE_NPZ)
    k = d["cluster_k"]
    sizes = np.array(sorted(set(k.tolist())))
    radius = np.array([d["cluster_radius"][k == kk].mean() for kk in sizes])

    rows = []
    for ion, key in (("EMIM+", "cluster_cation_q"), ("NTf2-", "cluster_anion_q")):
        q = np.abs(d[key])
        for kk, r in zip(sizes, radius):
            sample = q[k == kk]
            rows.append(dict(
                series="cluster", ion=ion, cluster_k=int(kk), radius_A=r,
                abs_q_e=sample.mean(),
                abs_q_sem_e=sample.std(ddof=1) / np.sqrt(sample.size),
                n_centres=sample.size))
    # The periodic-bulk value is a mean over all 250 ions of one frame, not over centres,
    # so it carries no SEM column: an error bar there would read as an uncertainty on the
    # reference line when it is really ion-to-ion variation.
    rows.append(dict(series="periodic_bulk", ion="EMIM+", cluster_k=pd.NA,
                     radius_A=pd.NA, abs_q_e=float(np.abs(d["bulk_cation_q"]).mean()),
                     abs_q_sem_e=pd.NA, n_centres=d["bulk_cation_q"].size))
    return pd.DataFrame(rows)


def neutral_tables(source: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Dipole ratios, and the geometries used for the structure thumbnails."""
    d = np.load(source / NEUTRAL_NPZ)
    rows = []
    for species, label in NEUTRALS:
        mu_dft, mu_so3lr = (float(x) for x in d[f"{species}_dipoles"])
        rows.append(dict(species=species, label=label, charged=False,
                         mu_dft_D=mu_dft, mu_so3lr_D=mu_so3lr,
                         ratio=mu_so3lr / mu_dft, n_geometries=1))

    # The ion-pair result averages the sampled 300 K geometries. Columns of
    # `SO3LR_geoms` are [q_MDCq, q_Mulliken, q_SO3LR, mu_DFT, mu_SO3LR]; the ratio is the
    # mean of the per-snapshot ratios, not the ratio of the means.
    pair = np.load(source / PAIR_NPZ)["SO3LR_geoms"]
    rows.append(dict(species=PAIR_SPECIES, label=PAIR_LABEL, charged=True,
                     mu_dft_D=float(pair[:, 3].mean()),
                     mu_so3lr_D=float(pair[:, 4].mean()),
                     ratio=float((pair[:, 4] / pair[:, 3]).mean()),
                     n_geometries=pair.shape[0]))

    geometries = []
    for species, _ in NEUTRALS:
        path = source / "ref_dft" / "PBE_TZP" / species / "ams.results" / "output.xyz"
        geometries.append((species, path))
    geometries.append((PAIR_SPECIES, source / PAIR_GEOMETRY))

    atoms = []
    for species, path in geometries:
        lines = path.read_text().splitlines()
        n = int(lines[0].split()[0])
        for index, line in enumerate(lines[2:2 + n]):
            element, x, y, z = line.split()[:4]
            atoms.append(dict(species=species, atom_index=index, element=element,
                              x_A=float(x), y_A=float(y), z_A=float(z)))
    return pd.DataFrame(rows), pd.DataFrame(atoms)


def pes_table(source: Path) -> pd.DataFrame:
    """The dissociation scan: relative energy and net cation charge per method.

    The conventional dE columns use separated fragments evaluated with the same method.
    SO3LR additionally has a pair-plateau-referenced curve, which does not use the
    separately evaluated formal-charge fragments.
    """
    c = np.load(source / PES_NPZ)
    dE_so3lr_pair_plateau = c["dE_ml_pair_plateau"]
    plateau_min_A = float(c["so3lr_pair_plateau_min_A"])
    table = pd.DataFrame({
        "separation_A": c["r"],
        "dE_pbe0_kJ_mol": c["dE_pbe0"],
        "dE_appleP_kJ_mol": c["dE_ap"],
        "dE_so3lr_kJ_mol": c["dE_ml"],
        "dE_so3lr_pair_plateau_kJ_mol": dE_so3lr_pair_plateau,
        "q_pbe0_mdcq_e": c["q_mdcq_pbe0"],
        "q_appleP_e": c["q_ap"],
        "q_so3lr_e": c["q_ml"],
    })
    table.attrs["so3lr_pair_plateau_min_A"] = plateau_min_A
    return table.sort_values("separation_A").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path, required=True,
        help="directory holding the compute_*.py output caches and DFT geometries")
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).resolve().parents[1] / "data",
                        help="chapter data/ directory to write the CSVs into")
    args = parser.parse_args()

    missing = [n for n in (CONVERGENCE_NPZ, NEUTRAL_NPZ, PAIR_NPZ, PES_NPZ)
               if not (args.source / n).is_file()]
    if missing:
        raise SystemExit(f"missing input caches in {args.source}: {', '.join(missing)}")

    args.out.mkdir(parents=True, exist_ok=True)

    convergence = convergence_table(args.source)
    path = args.out / "charge_convergence.csv"
    convergence.to_csv(path, index=False)
    _header(path, [
        "SO3LR net ionic charge versus charge-neutral cluster radius, 300 K bulk frame.",
        "series=cluster: mean |q| over 6 central ion pairs at each size k.",
        "series=periodic_bulk: mean |q| over all 250 cations of the periodic 8500-atom box.",
        "abs_q_sem_e is the SD over the centres divided by sqrt(n_centres); the bulk row",
        "has none because its scatter is ion-to-ion variation, not error on a mean.",
        f"The model's receptive field is {RECEPTIVE_FIELD_A} A = 3 message passes x 4.5 A.",
    ])
    print(f"wrote {path} ({len(convergence)} rows)")

    neutral, geometries = neutral_tables(args.source)
    path = args.out / "neutral_control.csv"
    neutral.to_csv(path, index=False)
    _header(path, [
        "Molecular dipole, SO3LR versus DFT at identical geometries.",
        "DFT: PBE-D3(BJ)/TZP, all-electron, geometry-optimised in AMS/ADF.",
        "The dipole is used because it needs no partitioning of the electron density.",
        "The ion-pair row uses 50 geometries.",
        "Its ratio is the mean of the per-geometry SO3LR/DFT ratios.",
    ])
    print(f"wrote {path} ({len(neutral)} rows)")

    path = args.out / "neutral_control_geometries.csv"
    geometries.to_csv(path, index=False)
    _header(path, [
        "Geometries behind the structure thumbnails, in Angstrom.",
        "Neutrals are DFT-optimised (PBE-D3/TZP); the ion-pair thumbnail is one",
        "representative of the 50 sampled geometries. Coordinates are as-computed; the notebook",
        "projects each molecule onto its own principal plane for display.",
    ])
    print(f"wrote {path} ({len(geometries)} atoms, "
          f"{geometries.species.nunique()} species)")

    pes = pes_table(args.source)
    plateau_min_A = pes.attrs["so3lr_pair_plateau_min_A"]
    path = args.out / "pes_scan.csv"
    pes.to_csv(path, index=False)
    _header(path, [
        "Gas-phase EMIM+/NTf2- dissociation scan at frozen internal geometry.",
        "separation_A is the centre-of-mass to centre-of-mass distance.",
        "The PBE0, APPLE&P, and dE_so3lr_kJ_mol columns use the sum of separately",
        "evaluated cation and anion energies from the same method as their zero.",
        "dE_so3lr_pair_plateau_kJ_mol instead uses the mean intact-pair energy for",
        f"separations >= {plateau_min_A:.1f} A, where the SO3LR scan is flat. This avoids",
        "using formal +/-1 isolated fragments for a distant pair whose SO3LR charges are",
        "not +/-1. The plateau-referenced curve is the physically relevant SO3LR curve.",
        "PBE0-D3/TZ2P is the electronic-structure reference.",
        "q columns are the net charge of EMIM+. APPLE&P's is fixed by parameterisation.",
    ])
    print(f"wrote {path} ({len(pes)} rows, "
          f"{pes.separation_A.min():.1f}-{pes.separation_A.max():.1f} A)")


if __name__ == "__main__":
    main()
