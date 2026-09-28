#!/usr/bin/env python3
"""Stream ion COM velocities and charge current from AMS RKF or SO3LR HDF5.

The external output contract is ``<prefix>_{ionvels,chargecurrent}.npy`` plus
metadata JSON. It never needs within-ion atom indices: every cation is a
complete 19-atom block and every anion a complete 15-atom block. Thus
APPLE&P's different raw order is handled safely by reading its own RKF masses.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
from numpy.lib.format import open_memmap

N_CAT = N_ANI = 250; CAT_ATOMS = 19; ANI_ATOMS = 15; N_ATOMS = 8500; ANI_START = 4750
JAX_VEL_TO_M_S = 9822.69641728187; BOHR_FS_TO_M_S = 5.29177210903e4
MASS = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998403, 'S': 32.06}
CAT_TEMPLATE = 'N C N C C C H C H H H H H C H H H H H'.split()
ANI_TEMPLATE = 'C F F F S N O O S O O C F F F'.split()


def close(tmp: Path, final: Path, array) -> None:
    array.flush()
    if getattr(array, '_mmap', None) is not None: array._mmap.close()
    os.replace(tmp, final)


def com(velocity: np.ndarray, masses: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    velocity = np.asarray(velocity, dtype=float).reshape(N_ATOMS, 3)
    cm = masses[:ANI_START].reshape(N_CAT, CAT_ATOMS); am = masses[ANI_START:].reshape(N_ANI, ANI_ATOMS)
    cation = np.einsum('iaj,ia->ij', velocity[:ANI_START].reshape(N_CAT, CAT_ATOMS, 3), cm / cm.sum(1)[:, None])
    anion = np.einsum('iaj,ia->ij', velocity[ANI_START:].reshape(N_ANI, ANI_ATOMS, 3), am / am.sum(1)[:, None])
    return cation, anion


def rkf_reader(path: Path):
    from scm.plams import RKFTrajectoryFile
    trajectory = RKFTrajectoryFile(str(path)); trajectory.store_mddata()
    section, item = trajectory.mdhistory_name, 'Velocities'
    if item in trajectory.mditems:
        read = lambda frame: trajectory.file_object.read(section, f'{item}({frame + 1})')
    elif item in trajectory.mdblockitems:
        cache = {}
        def read(frame):
            block, offset = frame // trajectory.mdblocksize + 1, frame % trajectory.mdblocksize
            if block not in cache:
                values = trajectory.file_object.read(section, f'{item}({block})')
                cache.clear(); cache[block] = values.split() if isinstance(values, str) else values
            return cache[block][offset]
    else:
        trajectory.close(); raise ValueError('RKF has no MDHistory Velocities')
    masses = np.array([atom.mass for atom in trajectory.get_plamsmol()], dtype=float)
    _, cell = trajectory.read_frame(0); cell = np.asarray(cell, dtype=float)
    box = np.diag(cell) if cell.shape == (3, 3) else cell
    return trajectory, read, trajectory.get_length(), masses, box, float(trajectory.timestep) / 1000.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', choices=('rkf', 'hdf5'), required=True)
    parser.add_argument('--force-field', choices=('apple', 'so3lr', 'clp', 'clp-scaled'), required=True)
    parser.add_argument('--rkf', type=Path); parser.add_argument('--hdf5', type=Path, nargs='+')
    parser.add_argument('--out-prefix', type=Path, required=True); parser.add_argument('--q-cation', type=float, default=1.0)
    parser.add_argument('--q-anion', type=float, default=-1.0); parser.add_argument('--saved-dt-ps', type=float, default=.005)
    parser.add_argument('--max-frames', type=int, default=0)
    args = parser.parse_args()
    if args.source == 'rkf' and not args.rkf: parser.error('--rkf is required')
    if args.source == 'hdf5' and not args.hdf5: parser.error('--hdf5 is required')
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    trajectory = None
    if args.source == 'rkf':
        trajectory, reader, n_frames, masses, box, dt_ps = rkf_reader(args.rkf); factor = BOHR_FS_TO_M_S; paths = [str(args.rkf)]
    else:
        import h5py
        masses = np.array([MASS[x] for x in CAT_TEMPLATE] * N_CAT + [MASS[x] for x in ANI_TEMPLATE] * N_ANI)
        n_frames = 0
        for number, path in enumerate(args.hdf5):
            with h5py.File(path, 'r') as handle: n_frames += handle['momenta'].shape[0] * handle['momenta'].shape[1] - bool(number)
        with h5py.File(args.hdf5[0], 'r') as handle: box = np.asarray(handle['box'][0, 0], dtype=float)
        dt_ps, factor, paths = args.saved_dt_ps, JAX_VEL_TO_M_S, [str(path) for path in args.hdf5]
    if args.max_frames: n_frames = min(n_frames, args.max_frames)
    ion_path = Path(f'{args.out_prefix}_ionvels.npy'); current_path = Path(f'{args.out_prefix}_chargecurrent.npy')
    ion_tmp, current_tmp = Path(str(ion_path) + '.tmp'), Path(str(current_path) + '.tmp')
    ionvels = open_memmap(ion_tmp, mode='w+', dtype=np.float32, shape=(n_frames, 500, 3))
    current = open_memmap(current_tmp, mode='w+', dtype=np.float64, shape=(n_frames, 3))
    def store(frame, velocity):
        cation, anion = com(velocity, masses); ionvels[frame, :250] = cation; ionvels[frame, 250:] = anion
        current[frame] = args.q_cation * cation.sum(0) + args.q_anion * anion.sum(0)
    if args.source == 'rkf':
        try:
            for frame in range(n_frames): store(frame, reader(frame))
        finally: trajectory.close()
    else:
        import h5py
        frame = 0
        for number, path in enumerate(args.hdf5):
            with h5py.File(path, 'r') as handle:
                for block in range(handle['momenta'].shape[0]):
                    momenta = np.asarray(handle['momenta'][block], dtype=float)
                    if number and block == 0: momenta = momenta[1:]
                    for momentum in momenta:
                        if frame >= n_frames: break
                        store(frame, momentum / masses[:, None]); frame += 1
                    if frame >= n_frames: break
            if frame >= n_frames: break
    close(ion_tmp, ion_path, ionvels); close(current_tmp, current_path, current)
    meta = {'force_field': args.force_field, 'source': args.source, 'source_paths': paths, 'n_frames': n_frames,
            'saved_dt_ps': dt_ps, 'volume_A3': float(np.prod(box)), 'q_cation_e': args.q_cation,
            'q_anion_e': args.q_anion, 'velocity_to_m_s': factor, 'ionvels': str(ion_path), 'chargecurrent': str(current_path)}
    Path(f'{args.out_prefix}_transport_metadata.json').write_text(json.dumps(meta, indent=2) + '\n')
    print(f'wrote {n_frames} frames to {ion_path} and {current_path}')


if __name__ == '__main__': main()
