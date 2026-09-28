#!/usr/bin/env python3
"""Create the compact transport plot table, including STACIE uncertainties."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def derived_errors(path: Path, force_field: str) -> pd.DataFrame:
    data = pd.read_csv(path).rename(columns={
        "T": "temperature_K", "D_cation_std": "D_cation_std_m2_s",
        "D_anion_std": "D_anion_std_m2_s", "sigma_std": "sigma_GK_std_S_m",
    })
    data["force_field"] = force_field
    return data[["force_field", "temperature_K", "D_cation_std_m2_s",
                 "D_anion_std_m2_s", "sigma_GK_std_S_m"]]


def stacie_errors(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path)
    rows = []
    for temperature, group in data.groupby("temperature_K"):
        record = {"force_field": "SO3LR", "temperature_K": temperature}
        for property_name, column in (("diffusion_cation", "D_cation_std_m2_s"),
                                      ("diffusion_anion", "D_anion_std_m2_s"),
                                      ("conductivity", "sigma_GK_std_S_m")):
            candidates = group[group.property.eq(property_name)]
            record[column] = candidates.loc[candidates.acint_std.idxmin(), "acint_std"]
        rows.append(record)
    return pd.DataFrame(rows)


def scaled_errors(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, comment="#").rename(columns={
        "T_K": "temperature_K", "D_cation_std_m2_s": "D_cation_std_m2_s",
        "D_anion_std_m2_s": "D_anion_std_m2_s",
        "sigma_std_S_per_m": "sigma_GK_std_S_m",
    })
    data["force_field"] = "CL&P-scaled"
    return data[["force_field", "temperature_K", "D_cation_std_m2_s",
                 "D_anion_std_m2_s", "sigma_GK_std_S_m"]]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--central", type=Path, required=True)
    parser.add_argument("--applep", type=Path, required=True)
    parser.add_argument("--so3lr", type=Path, required=True)
    parser.add_argument("--clp", type=Path, required=True)
    parser.add_argument("--clp-scaled", type=Path, required=True)
    parser.add_argument("--density", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    central = pd.read_csv(args.central)
    errors = pd.concat((derived_errors(args.applep, "APPLE&P"), stacie_errors(args.so3lr),
                        derived_errors(args.clp, "CL&P"), scaled_errors(args.clp_scaled)),
                       ignore_index=True)
    data = central.merge(errors, on=["force_field", "temperature_K"], validate="one_to_one")
    density = pd.read_csv(args.density).rename(columns={"density_sem_g_cm3": "density_std_g_cm3"})
    density["force_field"] = density.force_field.replace({"OPLS": "CL&P", "OPLS-scaled": "CL&P-scaled"})
    data = data.merge(density[["force_field", "temperature_K", "density_g_cm3", "density_std_g_cm3"]],
                      on=["force_field", "temperature_K"], validate="many_to_one")
    d_sum = data.D_cation_m2_s + data.D_anion_m2_s
    d_sum_std = np.hypot(data.D_cation_std_m2_s, data.D_anion_std_m2_s)
    data["sigma_NE_std_S_m"] = data.sigma_NE_S_m * np.hypot(d_sum_std / d_sum,
                                                               data.density_std_g_cm3 / data.density_g_cm3)
    data["haven_ratio_std"] = data.haven_ratio_GK_over_NE * np.hypot(
        data.sigma_GK_std_S_m / data.sigma_GK_S_m,
        data.sigma_NE_std_S_m / data.sigma_NE_S_m,
    )
    columns = ["force_field", "temperature_K", "D_cation_m2_s", "D_cation_std_m2_s",
               "D_anion_m2_s", "D_anion_std_m2_s", "sigma_GK_S_m", "sigma_GK_std_S_m",
               "sigma_NE_S_m", "sigma_NE_std_S_m", "haven_ratio_GK_over_NE", "haven_ratio_std"]
    data = data[columns].sort_values(["force_field", "temperature_K"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
