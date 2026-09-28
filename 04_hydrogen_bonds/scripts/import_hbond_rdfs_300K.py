#!/usr/bin/env python3
"""Export the compact 300 K all-ring-C--H RDF curves used by the figure."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def standard_curves(path: Path, force_field: str) -> list[pd.DataFrame]:
    result = np.load(path)
    r = result["r_A"]
    return [pd.DataFrame({"force_field": force_field, "pair": pair,
                          "r_A": r, "g_r": result[key]})
            for pair, key in (("H_O", "rdf_Hring_O"), ("H_F", "rdf_Hring_F"))]


def so3lr_curves(path: Path) -> list[pd.DataFrame]:
    result = np.load(path)
    return [pd.DataFrame({"force_field": "SO3LR", "pair": pair,
                          "r_A": result[f"SO3LR|ring|{atom}|r"],
                          "g_r": result[f"SO3LR|ring|{atom}|g"]})
            for pair, atom in (("H_O", "o"), ("H_F", "f"))]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--so3lr", type=Path, required=True,
                        help="matched_hbond_rdf.npz containing SO3LR ring curves")
    parser.add_argument("--applep", type=Path, required=True)
    parser.add_argument("--clp", type=Path, required=True)
    parser.add_argument("--clp-scaled", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    curves = so3lr_curves(args.so3lr)
    curves += standard_curves(args.applep, "APPLE&P")
    curves += standard_curves(args.clp, "CL&P")
    curves += standard_curves(args.clp_scaled, "CL&P-scaled")
    output = pd.concat(curves, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
