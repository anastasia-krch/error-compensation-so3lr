#!/usr/bin/env python3
"""Combine per-run STACIE JSON summaries and derive NE conductivity/Haven ratio."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.constants import Boltzmann, elementary_charge


N_IONS_PER_SPECIES = 250


def value(summary: dict, key: str) -> float:
    try:
        return float(summary[key]["value"])
    except KeyError as error:
        raise KeyError(f"missing {key}.value in STACIE summary") from error


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="per-run fit_stacie JSON files")
    parser.add_argument("--force-field", required=True, choices=("APPLE&P", "SO3LR", "CL&P", "CL&P-scaled"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    for path in args.inputs:
        summary = json.loads(path.read_text())
        temperature = float(summary["temperature_K"])
        volume_m3 = float(summary["volume_A3"]) * 1e-30
        q_e = abs(float(summary["cation_charge_e"]))
        d_cation = value(summary, "diffusion_cation_m2_s")
        d_anion = value(summary, "diffusion_anion_m2_s")
        sigma_gk = value(summary, "conductivity_S_m")
        sigma_ne = (N_IONS_PER_SPECIES * (q_e * elementary_charge) ** 2 *
                    (d_cation + d_anion) / (volume_m3 * Boltzmann * temperature))
        rows.append({"force_field": args.force_field, "temperature_K": temperature,
                     "D_cation_m2_s": d_cation, "D_anion_m2_s": d_anion,
                     "sigma_GK_S_m": sigma_gk, "sigma_NE_S_m": sigma_ne,
                     "haven_ratio_GK_over_NE": sigma_gk / sigma_ne,
                     "source_summary": str(path)})
    table = pd.DataFrame(rows).sort_values("temperature_K")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)
    print(f"wrote {args.output} ({len(table)} rows)")


if __name__ == "__main__":
    main()
