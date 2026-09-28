#!/usr/bin/env python3
"""Export all structural reduced trajectories from one EMIM-NTf2 trajectory pass.

Outputs (all asexyz with the instantaneous orthorhombic cell):
  * ``*_nmid_anionN.xyz``: 250 anion central N, then 250 cation ring-N
    midpoints.  Used for the N(midpoint)--N anion RDF and ion-pair work.
  * ``*_hbond.xyz``: ring C, their bonded ring H, NTf2 O, NTf2 F. Used by
    LionAnalysis H-bond RDF/lifetime configurations.
  * ``*_dihedral.xyz``: NTf2 C1, S1, N, S2, C2 atom groups. Used for the two
    backbone dihedrals C1-S1-N-S2 and S1-N-S2-C2.

The three force-field topology orders are explicit. APPLE&P is *not* compatible
with the SO3LR/CL&P order. Raw trajectories and reduced XYZ are deliberately
not publication-repository data; this script records enough provenance to
regenerate them in the simulation environment.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

N_CATION = N_ANION = 250
N_EMIM_ATOMS, N_NTF2_ATOMS = 19, 15
N_ATOMS = N_CATION * N_EMIM_ATOMS + N_ANION * N_NTF2_ATOMS
NTF2_START = N_CATION * N_EMIM_ATOMS

# Local offsets in the 19-atom EMIM and 15-atom NTf2 blocks, respectively.
# SO3LR and both CL&P variants use the same order; APPLE&P does not.
TOPOLOGY = {
    "so3lr": {
        "ring_ch": ((1, 6), (3, 8), (4, 9)), "ring_n": (0, 2),
        "anion_n": 5, "oxygen": (6, 7, 9, 10), "fluorine": (1, 2, 3, 12, 13, 14),
        "dihedral": (0, 4, 5, 8, 11),  # C1,S1,N,S2,C2
    },
    "clp": {
        "ring_ch": ((1, 6), (3, 8), (4, 9)), "ring_n": (0, 2),
        "anion_n": 5, "oxygen": (6, 7, 9, 10), "fluorine": (1, 2, 3, 12, 13, 14),
        "dihedral": (0, 4, 5, 8, 11),
    },
    "clp-scaled": {
        "ring_ch": ((1, 6), (3, 8), (4, 9)), "ring_n": (0, 2),
        "anion_n": 5, "oxygen": (6, 7, 9, 10), "fluorine": (1, 2, 3, 12, 13, 14),
        "dihedral": (0, 4, 5, 8, 11),
    },
    "apple": {
        "ring_ch": ((7, 18), (3, 13), (4, 14)), "ring_n": (2, 5),
        "anion_n": 3, "oxygen": (0, 2, 5, 6), "fluorine": (8, 9, 10, 12, 13, 14),
        "dihedral": (11, 1, 3, 4, 7),
    },
}


def atom_indices(force_field: str) -> dict[str, np.ndarray]:
    """Global atom indices for the selected topology; validates all group sizes."""
    top = TOPOLOGY[force_field]
    cbase = np.arange(N_CATION) * N_EMIM_ATOMS
    abase = NTF2_START + np.arange(N_ANION) * N_NTF2_ATOMS
    ring_c = np.concatenate([cbase + c for c, _ in top["ring_ch"]])
    ring_h = np.concatenate([cbase + h for _, h in top["ring_ch"]])
    result = {
        "ring_c": ring_c, "ring_h": ring_h,
        "oxygen": np.concatenate([abase + x for x in top["oxygen"]]),
        "fluorine": np.concatenate([abase + x for x in top["fluorine"]]),
        "ring_n1": cbase + top["ring_n"][0], "ring_n2": cbase + top["ring_n"][1],
        "anion_n": abase + top["anion_n"],
    }
    for label, offset in zip(("dih_c1", "dih_s1", "dih_n", "dih_s2", "dih_c2"), top["dihedral"]):
        result[label] = abase + offset
    expected = {"ring_c": 750, "ring_h": 750, "oxygen": 1000, "fluorine": 1500,
                "ring_n1": 250, "ring_n2": 250, "anion_n": 250,
                "dih_c1": 250, "dih_s1": 250, "dih_n": 250, "dih_s2": 250, "dih_c2": 250}
    for label, size in expected.items():
        assert len(result[label]) == size and np.all((result[label] >= 0) & (result[label] < N_ATOMS))
    return {key: value.astype(np.int64) for key, value in result.items()}


def minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    return delta - box * np.rint(delta / box)


def ring_n_midpoint(positions: np.ndarray, box: np.ndarray, n1: np.ndarray, n2: np.ndarray) -> np.ndarray:
    first = positions[n1]
    return (first + 0.5 * minimum_image(positions[n2] - first, box)) % box


def write_frame(handle, groups, box: np.ndarray, time_ps: float, frame: int) -> None:
    """Write a reduced asexyz frame; group order is the documented ID order."""
    count = sum(len(coords) for _, coords in groups)
    lattice = f'{box[0]:.8f} 0 0 0 {box[1]:.8f} 0 0 0 {box[2]:.8f}'
    handle.write(f"{count}\nLattice=\"{lattice}\" Properties=species:S:1:pos:R:3 Time={time_ps:.6f} saved_frame={frame} pbc=\"T T T\"\n")
    for element, coordinates in groups:
        for x, y, z in coordinates:
            handle.write(f"{element} {x:.8f} {y:.8f} {z:.8f}\n")


def wanted_indices(total: int, saved_dt: float, hb_dt: float, hb_window: float,
                   structure_dt: float) -> set[int]:
    def stride(step: float) -> int:
        value = round(step / saved_dt)
        if value < 1 or not np.isclose(value * saved_dt, step):
            raise ValueError(f"{step} ps is not an integer multiple of saved dt {saved_dt} ps")
        return int(value)
    hb_stride, structure_stride = stride(hb_dt), stride(structure_dt)
    hb_end = total if hb_window <= 0 else min(total, int(round(hb_window / saved_dt)))
    return set(range(0, hb_end, hb_stride)) | set(range(0, total, structure_stride))


def rkf_frames(path: Path, indices: set[int]):
    from scm.plams import RKFTrajectoryFile
    trajectory = RKFTrajectoryFile(str(path))
    try:
        _, first_cell = trajectory.read_frame(0)
        fallback_box = np.diag(np.asarray(first_cell, dtype=float))
        for index in sorted(indices):
            positions, cell = trajectory.read_frame(index)
            box = fallback_box if cell is None else np.diag(np.asarray(cell, dtype=float))
            yield index, np.asarray(positions, dtype=float), box
    finally:
        trajectory.close()


def hdf5_frames(paths: list[Path], indices: set[int]):
    import h5py
    target = set(indices)
    frame = 0
    for path in paths:
        with h5py.File(path, "r") as handle:
            positions = handle["positions"]
            boxes = handle["box"]
            for chunk in range(positions.shape[0]):
                xyz_chunk = np.asarray(positions[chunk], dtype=float)
                box_chunk = np.asarray(boxes[chunk], dtype=float)
                for local in range(xyz_chunk.shape[0]):
                    if frame in target:
                        box = box_chunk[local]
                        box = np.diag(box) if box.shape == (3, 3) else box
                        xyz = xyz_chunk[local]
                        # The SO3LR files may store fractional coordinates.
                        if np.nanmax(xyz) <= 1.5:
                            xyz = xyz * box
                        yield frame, xyz, box
                    frame += 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=("rkf", "hdf5"), required=False)
    parser.add_argument("--rkf", type=Path, help="AMS RKF trajectory for APPLE&P/CL&P")
    parser.add_argument("--hdf5", type=Path, nargs="+", help="SO3LR HDF5 part(s)")
    parser.add_argument("--force-field", choices=tuple(TOPOLOGY), required=False)
    parser.add_argument("--out-prefix", type=Path)
    parser.add_argument("--saved-dt-ps", type=float, default=0.005)
    parser.add_argument("--hbond-dt-ps", type=float, default=0.05)
    parser.add_argument("--hbond-window-ps", type=float, default=100.0)
    parser.add_argument("--structure-dt-ps", type=float, default=1.0,
                        help="spacing for RDF/ion-pair and dihedral XYZs")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        for name in TOPOLOGY:
            atom_indices(name)
        print("topology self-test passed for APPLE&P, SO3LR, CL&P, and CL&P-scaled")
        return
    if not (args.source and args.force_field and args.out_prefix):
        parser.error("--source, --force-field, and --out-prefix are required")
    if args.source == "rkf" and not args.rkf:
        parser.error("--rkf is required with --source rkf")
    if args.source == "hdf5" and not args.hdf5:
        parser.error("--hdf5 is required with --source hdf5")

    if args.source == "rkf":
        from scm.plams import RKFTrajectoryFile
        probe = RKFTrajectoryFile(str(args.rkf)); total = probe.get_length(); probe.close()
    else:
        import h5py
        total = sum(h5py.File(path, "r")["positions"].shape[0] * h5py.File(path, "r")["positions"].shape[1] for path in args.hdf5)
    selected = wanted_indices(total, args.saved_dt_ps, args.hbond_dt_ps, args.hbond_window_ps, args.structure_dt_ps)
    idx = atom_indices(args.force_field)
    prefix = args.out_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    paths = {"nmid_anionN": Path(f"{prefix}_nmid_anionN.xyz"), "hbond": Path(f"{prefix}_hbond.xyz"),
             "dihedral": Path(f"{prefix}_dihedral.xyz")}
    hb_stride = round(args.hbond_dt_ps / args.saved_dt_ps)
    structure_stride = round(args.structure_dt_ps / args.saved_dt_ps)
    hb_end = total if args.hbond_window_ps <= 0 else min(total, round(args.hbond_window_ps / args.saved_dt_ps))
    frames = rkf_frames(args.rkf, selected) if args.source == "rkf" else hdf5_frames(args.hdf5, selected)
    counts = {key: 0 for key in paths}
    with paths["nmid_anionN"].open("w") as nmid, paths["hbond"].open("w") as hbond, paths["dihedral"].open("w") as dih:
        for frame, positions, box in frames:
            if positions.shape != (N_ATOMS, 3):
                raise ValueError(f"frame {frame}: expected ({N_ATOMS},3), got {positions.shape}")
            wrapped = positions % box
            time_ps = frame * args.saved_dt_ps
            if frame % structure_stride == 0:
                midpoint = ring_n_midpoint(wrapped, box, idx["ring_n1"], idx["ring_n2"])
                write_frame(nmid, [("N", wrapped[idx["anion_n"]]), ("N", midpoint)], box, time_ps, frame)
                write_frame(dih, [("C", wrapped[idx["dih_c1"]]), ("S", wrapped[idx["dih_s1"]]), ("N", wrapped[idx["dih_n"]]),
                                  ("S", wrapped[idx["dih_s2"]]), ("C", wrapped[idx["dih_c2"]])], box, time_ps, frame)
                counts["nmid_anionN"] += 1; counts["dihedral"] += 1
            if frame < hb_end and frame % hb_stride == 0:
                write_frame(hbond, [("C", wrapped[idx["ring_c"]]), ("H", wrapped[idx["ring_h"]]),
                                    ("O", wrapped[idx["oxygen"]]), ("F", wrapped[idx["fluorine"]])], box, time_ps, frame)
                counts["hbond"] += 1
    summary = {"force_field": args.force_field, "source": args.source, "saved_dt_ps": args.saved_dt_ps,
               "structure_dt_ps": args.structure_dt_ps, "hbond_dt_ps": args.hbond_dt_ps,
               "hbond_window_ps": args.hbond_window_ps, "n_frames": counts,
               "layouts_1based": {"nmid_anionN": "anion N 1-250; cation ring-N midpoint 251-500",
                 "hbond": "ring C 1-750; ring H 751-1500; NTf2 O 1501-2500; NTf2 F 2501-4000",
                 "dihedral": "C1 1-250; S1 251-500; N 501-750; S2 751-1000; C2 1001-1250"}}
    Path(f"{prefix}_export_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"wrote {counts} reduced-trajectory frames")


if __name__ == "__main__":
    main()
