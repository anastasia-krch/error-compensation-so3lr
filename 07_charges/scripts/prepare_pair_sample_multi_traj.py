#!/usr/bin/env python3
"""Sample EMIM+/NTf2- ion-pair geometries from independent NVE trajectories.

Ten pairs are selected from one frame of each trajectory by default. With five input
trajectories, this produces the 50 geometries used for the dipole comparison.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_charge_convergence import build_cluster, choose_centres  # noqa: E402

# Atom order recovered and cross-checked in get_vel_hdf5.py (the HDF5 stores no species):
# 250 EMIM+ (19 atoms each) followed by 250 NTf2- (15 atoms each).
EMIM = ["N", "C", "N", "C", "C", "C", "H", "C", "H", "H",
        "H", "H", "H", "C", "H", "H", "H", "H", "H"]
NTF2 = ["C", "F", "F", "F", "S", "N", "O", "O", "S", "O",
        "O", "C", "F", "F", "F"]
N_EMIM = N_NTF2 = 250


def read_frame(hdf5_path: Path, frame_index: int):
    """One (positions, box) pair, touching only the one HDF5 chunk that holds it.

    `positions` is stored as fractional/scaled coordinates in [0, 1), not Angstrom -- a
    frame read directly and treated as Cartesian puts every atom within ~1 Angstrom of
    every other, which turns fragment detection into a near-complete bond graph and hangs.
    Scale by the (orthorhombic, diagonal) box to get Cartesian Angstrom.
    """
    import h5py

    with h5py.File(hdf5_path, "r") as f:
        per_block = f["positions"].shape[1]
        block, offset = divmod(frame_index, per_block)
        fractional = np.asarray(f["positions"][block, offset], dtype=np.float64)
        box = np.asarray(f["box"][block, offset], dtype=np.float64)
    if not np.isfinite(fractional).all():
        raise RuntimeError(f"{hdf5_path} frame {frame_index}: non-finite positions")
    positions = fractional * box[None, :]
    return positions, box


def build_atoms(positions, box):
    from ase import Atoms

    symbols = EMIM * N_EMIM + NTF2 * N_NTF2
    return Atoms(symbols=symbols, positions=positions, cell=np.diag(box), pbc=True)


def fixed_fragments():
    """Cation/anion index groups from the HDF5's known, fixed atom ordering.

    Unlike charge_common.detect_ion_fragments (bond-distance based, for AMS output where
    atom identity isn't known in advance), this trajectory's atom order is deterministic:
    250 cations of 19 atoms then 250 anions of 15 atoms, per get_vel_hdf5.py's mass array.
    Distance-based bonding is best avoided here anyway -- in a 250-anion box, some
    sulfonyl...sulfonyl contact between two UNRELATED anions is close enough to trip a
    natural_cutoffs bond on many frames, which has nothing to do with the pair being
    sampled but fails a global "exactly two fragment sizes" check.
    """
    cations = [np.arange(m * len(EMIM), (m + 1) * len(EMIM)) for m in range(N_EMIM)]
    offset = N_EMIM * len(EMIM)
    anions = [np.arange(offset + m * len(NTF2), offset + (m + 1) * len(NTF2))
              for m in range(N_NTF2)]
    return anions, cations


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hdf5", type=Path, nargs="+", required=True,
                        help="independent trajectory files, one pair-set drawn per file")
    parser.add_argument("--frame-fraction", type=float, default=0.5,
                        help="fraction into each trajectory to draw the shared frame from")
    parser.add_argument("--n-centres", type=int, default=10,
                        help="ion pairs sampled per trajectory")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    manifest_rows = []

    for traj_id, hdf5_path in enumerate(args.hdf5):
        import h5py

        with h5py.File(hdf5_path, "r") as f:
            n_blocks, per_block = f["positions"].shape[:2]
        frame_index = int(round(args.frame_fraction * (n_blocks * per_block - 1)))
        positions, box = read_frame(hdf5_path, frame_index)
        atoms = build_atoms(positions, box)
        anions, cations = fixed_fragments()
        print(f"{hdf5_path.parent.name}: frame {frame_index}/{n_blocks * per_block}, "
              f"{len(cations)} cations, {len(anions)} anions, box {box}")

        for cation, anion in choose_centres(atoms, anions, cations, args.n_centres, seed=traj_id):
            cluster, _radius = build_cluster(atoms, anions, cations, (cation, anion), k=0)
            tag = f"pair_{hdf5_path.parent.name}_c{cation:03d}_a{anion:03d}"
            path = args.out / f"{tag}.xyz"
            cluster.write(path)
            manifest_rows.append((f"{tag}.xyz", 0, tag))

    manifest = args.out / "pair_manifest_multi_traj.tsv"
    with manifest.open("w") as handle:
        handle.write("xyz\tcharge\ttag\n")
        for filename, charge, tag in manifest_rows:
            handle.write(f"{filename}\t{charge}\t{tag}\n")
    print(f"wrote {manifest} ({len(manifest_rows)} pairs from {len(args.hdf5)} trajectories)")


if __name__ == "__main__":
    main()
