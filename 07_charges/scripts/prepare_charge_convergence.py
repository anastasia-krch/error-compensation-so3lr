#!/usr/bin/env python3
"""Prepare frozen charge-neutral clusters for manifest-driven AMS/SO3LR runs.

This is deliberately separate from ``compute_charge_convergence.py``: this script only
selects and writes geometries, whereas AMS evaluates every SO3LR charge. The two stages
make the cluster analysis reproducible through the same AMS driver used for DFT and
APPLE&P, and let a Slurm array retry one failed structure without repeating JAX work.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from charge_common import N_LAYERS, SR_CUTOFF_A, detect_ion_fragments, unwrap_fragment


DEFAULT_K_VALUES = (0, 1, 2, 3, 4, 6, 8, 12, 16, 20, 26, 32, 48, 64, 96)


def build_cluster(atoms, anions, cations, centre_pair, k):
    """Central pair plus the k nearest cations and k nearest anions, at fixed geometry."""
    from ase import Atoms

    cell_diagonal = np.diag(np.asarray(atoms.get_cell()))
    if not np.all(cell_diagonal > 0):
        raise RuntimeError("charge-convergence input needs an orthorhombic periodic cell")
    positions = atoms.get_positions()
    centre_cation, centre_anion = centre_pair

    whole = {}
    for label, fragments in (("cation", cations), ("anion", anions)):
        for index, fragment in enumerate(fragments):
            whole[(label, index)] = unwrap_fragment(positions, fragment, cell_diagonal)

    pair_cation = whole[("cation", centre_cation)]
    pair_anion = whole[("anion", centre_anion)]
    pair_anion = pair_anion - np.round(
        (pair_anion.mean(axis=0) - pair_cation.mean(axis=0)) / cell_diagonal) * cell_diagonal
    origin = np.vstack((pair_cation, pair_anion)).mean(axis=0)

    def nearest_image(key):
        coordinates = whole[key]
        return coordinates - np.round((coordinates.mean(axis=0) - origin) / cell_diagonal) * cell_diagonal

    def ranked(label, fragments, skip):
        distances = [(np.linalg.norm(nearest_image((label, index)).mean(axis=0) - origin), index)
                     for index in range(len(fragments)) if index != skip]
        return [index for _, index in sorted(distances)]

    keys = [("cation", centre_cation), ("anion", centre_anion)]
    keys.extend(("cation", index) for index in ranked("cation", cations, centre_cation)[:k])
    keys.extend(("anion", index) for index in ranked("anion", anions, centre_anion)[:k])

    numbers, coordinates = [], []
    atomic_numbers = atoms.get_atomic_numbers()
    for label, index in keys:
        fragment = (cations if label == "cation" else anions)[index]
        numbers.append(atomic_numbers[fragment])
        coordinates.append(nearest_image((label, index)))
    cluster = Atoms(numbers=np.concatenate(numbers), positions=np.vstack(coordinates), pbc=False)
    radius = np.linalg.norm(cluster.get_positions() - origin, axis=1).max()
    return cluster, radius


def choose_centres(atoms, anions, cations, n_centres, seed):
    """Sample representative cations and pair each with its nearest anion."""
    cell_diagonal = np.diag(np.asarray(atoms.get_cell()))
    positions = atoms.get_positions()
    cation_com = np.array([unwrap_fragment(positions, fragment, cell_diagonal).mean(axis=0)
                           for fragment in cations])
    anion_com = np.array([unwrap_fragment(positions, fragment, cell_diagonal).mean(axis=0)
                          for fragment in anions])
    delta = anion_com[None, :, :] - cation_com[:, None, :]
    delta -= np.round(delta / cell_diagonal) * cell_diagonal
    nearest = np.linalg.norm(delta, axis=-1).argmin(axis=1)
    chosen = np.random.default_rng(seed).choice(len(cation_com), size=n_centres, replace=False)
    return [(int(cation), int(nearest[cation])) for cation in chosen]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frame", type=Path, required=True,
                        help="equilibrated periodic frame (extended XYZ with a lattice)")
    parser.add_argument("--out", type=Path, required=True,
                        help="directory for the frozen XYZs, manifest and metadata")
    parser.add_argument("--n-centres", type=int, default=6)
    parser.add_argument("--k-values", type=int, nargs="+", default=list(DEFAULT_K_VALUES))
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--skip-bulk", action="store_true",
                        help="do not write the full periodic reference geometry")
    args = parser.parse_args()

    from ase.io import read, write

    atoms = read(args.frame)
    anions, cations = detect_ion_fragments(atoms)
    args.out.mkdir(parents=True, exist_ok=True)
    print(f"{len(atoms)} atoms, {len(cations)} cations; receptive field "
          f"{SR_CUTOFF_A * N_LAYERS:.1f} A")

    tags, k_values, radii, filenames = [], [], [], []
    for cation, anion in choose_centres(atoms, anions, cations, args.n_centres, args.seed):
        for k in args.k_values:
            cluster, radius = build_cluster(atoms, anions, cations, (cation, anion), k)
            tag = f"cluster_c{cation:03d}_a{anion:03d}_k{k:03d}"
            filename = f"{tag}.xyz"
            write(args.out / filename, cluster)
            tags.append(tag)
            k_values.append(k)
            radii.append(radius)
            filenames.append(filename)

    if not args.skip_bulk:
        write(args.out / "bulk_300K.xyz", atoms)
        tags.append("bulk_300K")
        k_values.append(-1)
        radii.append(np.nan)
        filenames.append("bulk_300K.xyz")

    manifest = args.out / "cluster_manifest.tsv"
    with manifest.open("w") as handle:
        handle.write("xyz\tcharge\ttag\n")
        for filename, tag in zip(filenames, tags):
            handle.write(f"{filename}\t0\t{tag}\n")
    metadata = args.out / "cluster_metadata.npz"
    np.savez(metadata, tag=np.array(tags), cluster_k=np.array(k_values),
             cluster_radius=np.array(radii), xyz=np.array(filenames),
             bulk_cation_indices=np.asarray(cations, dtype=int))
    print(f"wrote {manifest} and {metadata} ({len(tags)} AMS jobs)")


if __name__ == "__main__":
    main()
