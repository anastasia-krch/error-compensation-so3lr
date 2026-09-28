#!/usr/bin/env python3
"""Collect cluster-convergence charges from AMS/SO3LR single points.

Run ``prepare_charge_convergence.py`` first, then use its ``cluster_manifest.tsv`` with
``run_ams_singlepoint.sh`` and ``ENGINE=so3lr``. No direct SO3LR/JAX calculation occurs
here: every value is read from the AMS log produced by the charge-aware ASE adapter.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from charge_common import ANION_SIZE, CATION_SIZE, index_ams_logs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True,
                        help="cluster_metadata.npz written by prepare_charge_convergence.py")
    parser.add_argument("--logdir", type=Path, required=True,
                        help="AMS stdout logs from ENGINE=so3lr")
    parser.add_argument("--level", default="SO3LR/ASE")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    metadata = np.load(args.metadata)
    index = index_ams_logs(args.logdir)
    if not index:
        raise SystemExit(f"no AMS logs with provenance lines found in {args.logdir}")

    rows, bulk, missing = [], None, []
    for tag, k, radius in zip(metadata["tag"], metadata["cluster_k"], metadata["cluster_radius"]):
        tag = str(tag)
        result = index.get((tag, args.level))
        if result is None or result["atomic_charges"] is None:
            missing.append(tag)
            continue
        charges = np.asarray(result["atomic_charges"], dtype=float)
        if not np.isfinite(charges).all() or abs(charges.sum()) > 1e-5:
            raise SystemExit(f"{tag}: invalid SO3LR charges (sum={charges.sum():+.3e})")
        if tag == "bulk_300K":
            cation_indices = np.asarray(metadata["bulk_cation_indices"], dtype=int)
            expected = int(cation_indices.max()) + 1
            if len(charges) != expected:
                raise SystemExit(f"{tag}: expected {expected} atomic charges, found {len(charges)}")
            bulk = charges[cation_indices].sum(axis=1)
            continue
        if len(charges) < CATION_SIZE + ANION_SIZE:
            raise SystemExit(f"{tag}: cluster has only {len(charges)} atomic charges")
        rows.append((int(k), float(radius), charges[:CATION_SIZE].sum(),
                     charges[CATION_SIZE:CATION_SIZE + ANION_SIZE].sum()))
    if missing:
        raise SystemExit("missing usable SO3LR AMS results for: " + ", ".join(missing))
    if not rows:
        raise SystemExit("no cluster rows found in metadata")
    if bulk is None:
        bulk = np.full(len(metadata["bulk_cation_indices"]), np.nan)

    rows.sort(key=lambda row: (row[0], row[1]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output,
             cluster_k=np.array([row[0] for row in rows]),
             cluster_radius=np.array([row[1] for row in rows]),
             cluster_cation_q=np.array([row[2] for row in rows]),
             cluster_anion_q=np.array([row[3] for row in rows]),
             bulk_cation_q=bulk)
    print(f"wrote {args.output} ({len(rows)} cluster evaluations; {len(bulk)} bulk cations)")


if __name__ == "__main__":
    main()
