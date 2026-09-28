#!/usr/bin/env python3
"""Plot direct first-1/e nearest-neighbour lifetimes from a seed-level CSV.

Input columns are force_field, temperature_K, seed, tau_1_over_e_ps,
max_time_ps, and C_at_max_time. Blank tau values are right-censored runs. A
temperature mean and its sample standard deviation are reported only when all
included seeds crossed 1/e. No parametric fit is performed.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

STYLE = {"apple": ("APPLE&P", "#0077BB", "o"), "so3lr": ("SO3LR", "#EE7733", "s"),
         "opls": ("CL&P", "#009988", "^"), "scaled_opls": ("CL&P-scaled", "#CC3311", "D")}


def number(value):
    return float(value) if value not in (None, "", "nan", "NaN") else float("nan")


def load(path):
    with path.open(newline="") as handle:
        rows = [{**row, "temperature_K": int(row["temperature_K"]),
                 "seed": int(row["seed"]), "tau_1_over_e_ps": number(row["tau_1_over_e_ps"]),
                 "max_time_ps": number(row["max_time_ps"])} for row in csv.DictReader(handle)]
    if not rows:
        raise ValueError(f"{path}: no seed-level rows")
    return rows


def summarize(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["force_field"], row["temperature_K"])].append(row)
    result = []
    for (ff, temp), group in sorted(grouped.items()):
        finite = [row["tau_1_over_e_ps"] for row in group if np.isfinite(row["tau_1_over_e_ps"])]
        censored = [row["max_time_ps"] for row in group if not np.isfinite(row["tau_1_over_e_ps"])]
        complete = len(finite) == len(group)
        result.append({"force_field": ff, "temperature_K": temp, "n_seeds": len(group),
                       "n_crossed": len(finite),
                       "tau_1_over_e_mean_ps": np.mean(finite) if complete else np.nan,
                       "tau_1_over_e_seed_sd_ps": np.std(finite, ddof=1) if complete and len(finite) > 1 else np.nan,
                       "tau_1_over_e_lower_bound_ps": min(censored) if censored else np.nan})
    return result


def plot(rows, out, logy):
    plt.rcParams.update({'font.size': 14, 'font.family': 'sans-serif', 'axes.labelsize': 16,
                         'xtick.labelsize': 14, 'ytick.labelsize': 14, 'axes.linewidth': 1.5,
                         'lines.linewidth': 2.5, 'lines.markersize': 8})
    fig, ax = plt.subplots(figsize=(9, 7))
    for ff, (label, color, marker) in STYLE.items():
        selected = [row for row in rows if row["force_field"] == ff]
        crossed = [row for row in selected if np.isfinite(row["tau_1_over_e_mean_ps"])]
        if crossed:
            ax.errorbar([r["temperature_K"] for r in crossed], [r["tau_1_over_e_mean_ps"] for r in crossed],
                        yerr=[0 if not np.isfinite(r["tau_1_over_e_seed_sd_ps"]) else r["tau_1_over_e_seed_sd_ps"] for r in crossed],
                        color=color, marker=marker, linestyle="-", capsize=5, capthick=1.5,
                        elinewidth=1.4, label=label)
        for row in selected:
            if np.isfinite(row["tau_1_over_e_lower_bound_ps"]):
                y = row["tau_1_over_e_lower_bound_ps"]
                ax.plot(row["temperature_K"], y, marker=marker, color=color, ms=8, linestyle="none")
                ax.annotate(f">{y:.0f}", (row["temperature_K"], y), xytext=(5, 6), textcoords="offset points", fontsize=11)
    ax.set(xlim=(288, 333), xticks=(290, 300, 310, 320, 330), xlabel="Temperature (K)",
           ylabel=r"Ion-pair lifetime, $\tau_{1/e}$ (ps)")
    if logy:
        ax.set_yscale("log")
    for spine in ax.spines.values(): spine.set_linewidth(1.5)
    ax.tick_params(width=1.5, length=7); ax.legend(frameon=False, fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = summarize(load(args.seed_csv))
    with (args.out_dir / "nearest_neighbor_ionpair_tau_1e_all_forcefields.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0])); writer.writeheader(); writer.writerows(summary)
    plot(summary, args.out_dir / "nearest_neighbor_ionpair_tau_1e_all_forcefields_logy.png", True)


if __name__ == "__main__":
    main()
