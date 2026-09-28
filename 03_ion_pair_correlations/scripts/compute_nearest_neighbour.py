#!/usr/bin/env python3
"""Cation-centred, re-forming nearest-neighbour ion-pair correlation.

The reduced extended-XYZ input contains 250 anion central-N sites followed by
250 cation ring-N-midpoint sites in each periodic orthorhombic frame. For each
cation and frame, the nearest anion identity is assigned. The intermittent
correlation permits re-formation and is reported through its direct first
C(t)=exp(-1) crossing; no fit or distance cutoff is used.
"""
from __future__ import annotations

import argparse
import itertools
import json
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

N_ANION = N_CATION = 250
N_ATOMS = N_ANION + N_CATION


class IntermittentCorrelator:
    def __init__(self, max_lag):
        self.max_lag = max_lag
        self.inter = np.zeros(max_lag + 1)
        self.counts = []
        self.history = defaultdict(deque)

    def update(self, active, sample):
        self.counts.append(len(active))
        self.inter[0] += len(active)
        for pair in active:
            q = self.history[pair]
            while q and q[0] < sample - self.max_lag:
                q.popleft()
            for origin in q:
                self.inter[sample - origin] += 1
            q.append(sample)

    def finish(self):
        counts = np.asarray(self.counts, dtype=float)
        prefix = np.r_[0.0, np.cumsum(counts)]
        denom = np.array([prefix[len(counts) - lag] if len(counts) > lag else 0.0
                          for lag in range(self.max_lag + 1)])
        return np.divide(
            self.inter, denom, out=np.full_like(self.inter, np.nan), where=denom > 0)


def frames(path):
    with path.open() as handle:
        while True:
            first = handle.readline()
            if not first:
                return
            if int(first) != N_ATOMS:
                raise ValueError(f"{path}: expected {N_ATOMS} sites per frame")
            header = handle.readline()
            try:
                lattice = np.fromstring(header.split('Lattice="')[1].split('"')[0], sep=" ")
            except IndexError as exc:
                raise ValueError(f"{path}: frame lacks an extended-XYZ Lattice field") from exc
            if lattice.size != 9:
                raise ValueError(f"{path}: expected a 3x3 cell")
            box = lattice[[0, 4, 8]]
            if np.any(box <= 0):
                raise ValueError(f"{path}: invalid orthorhombic cell {box}")
            pos = np.empty((N_ATOMS, 3))
            for index in range(N_ATOMS):
                line = handle.readline().split()
                if len(line) < 4:
                    raise ValueError(f"{path}: truncated coordinate frame")
                pos[index] = [float(value) for value in line[1:4]]
            yield pos[:N_ANION], pos[N_ANION:], box


def first_crossing(time, corr):
    below = np.flatnonzero(corr <= np.exp(-1.0))
    if not len(below):
        return None
    i = int(below[0])
    if i == 0 or corr[i] == corr[i - 1]:
        return float(time[i])
    return float(time[i - 1] + (np.exp(-1.0) - corr[i - 1]) *
                 (time[i] - time[i - 1]) / (corr[i] - corr[i - 1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xyz", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--dt-ps", type=float, default=1.0)
    parser.add_argument("--max-time-ps", type=float,
                        help="truncate to this inclusive common-window time")
    args = parser.parse_args()
    if args.dt_ps <= 0:
        parser.error("--dt-ps must be positive")

    if args.max_time_ps is not None and args.max_time_ps < 0:
        parser.error("--max-time-ps must be non-negative")

    source = frames(args.xyz)
    if args.max_time_ps is not None:
        n_frames = int(np.floor(args.max_time_ps / args.dt_ps + 1e-12)) + 1
        source = itertools.islice(source, n_frames)
    saved = list(source)
    if len(saved) < 2:
        raise SystemExit("at least two frames are required")
    acc = IntermittentCorrelator(len(saved) - 1)
    for sample, (anion, cation, box) in enumerate(saved):
        tree = cKDTree(anion % box, boxsize=box)
        _, nearest = tree.query(cation % box, k=1)
        acc.update({i * N_ANION + int(j) for i, j in enumerate(nearest)}, sample)
    intermittent = acc.finish()
    time = np.arange(len(intermittent), dtype=float) * args.dt_ps
    args.out_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.out_dir / f"nn_pair_lifetime_intermittent_{args.label}.txt",
               np.c_[time, intermittent], header="time_ps C_intermittent")
    tau = first_crossing(time, intermittent)
    summary = {"label": args.label, "n_frames": len(saved), "dt_ps": args.dt_ps,
               "tau_1_over_e_ps": tau, "right_censored": tau is None,
               "max_time_ps": float(time[-1]), "C_at_max_time": float(intermittent[-1])}
    (args.out_dir / f"nn_pair_lifetime_summary_{args.label}.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
