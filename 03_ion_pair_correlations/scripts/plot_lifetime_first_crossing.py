#!/usr/bin/env python3
"""Render the paper example of the nearest-neighbour lifetime definition."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


CHAPTER = Path(__file__).resolve().parents[1]
SOURCES = sorted((CHAPTER / "data").glob("example_Ct_clp_scaled_310K_seed*.txt"))
OUTPUT = CHAPTER / "figures" / "nearest_neighbor_lifetime_first_crossing_example.png"
TARGET = np.exp(-1.0)


def first_1e_crossing(time_ps, correlation):
    """Linearly interpolate the first downward C(t)=exp(-1) crossing."""
    hit = np.flatnonzero(correlation <= TARGET)
    if not hit.size:
        raise ValueError("An example correlation does not reach exp(-1).")
    index = int(hit[0])
    if index == 0:
        return float(time_ps[0])
    t0, t1 = time_ps[index - 1], time_ps[index]
    c0, c1 = correlation[index - 1], correlation[index]
    return float(t0 + (TARGET - c0) * (t1 - t0) / (c1 - c0))


def main():
    if len(SOURCES) != 5:
        raise ValueError(f"Expected five repeat curves, found {len(SOURCES)}.")
    loaded = [np.loadtxt(path, comments="#") for path in SOURCES]
    time_ps = loaded[0][:, 0]
    correlations = np.array([curve[:, 1] for curve in loaded])
    if any(not np.array_equal(curve[:, 0], time_ps) for curve in loaded[1:]):
        raise ValueError("Repeat curves do not share the same time grid.")
    taus = np.array([first_1e_crossing(time_ps, curve) for curve in correlations])
    mean_tau = float(np.mean(taus))
    sd_tau = float(np.std(taus, ddof=1))
    xmax = 25 * np.ceil(1.5 * np.max(taus) / 25)
    colour = "#c83817"  # CL&P-scaled, shared density/transport palette

    plt.rcParams.update({
        "font.size": 14, "font.family": "sans-serif", "axes.labelsize": 16,
        "xtick.labelsize": 14, "ytick.labelsize": 14, "axes.linewidth": 1.5,
        "lines.linewidth": 2.5,
    })
    fig, ax = plt.subplots(figsize=(9, 6.5))
    shown = time_ps <= xmax
    for correlation in correlations:
        ax.plot(time_ps[shown], correlation[shown], color=colour, alpha=.30, lw=1.6)
    ax.plot(time_ps[shown], np.mean(correlations, axis=0)[shown], color=colour, lw=2.8)
    ax.axhline(TARGET, color="0.35", ls="--", lw=1.5, zorder=0)
    ax.scatter(taus, np.full_like(taus, TARGET), color=colour, s=42, zorder=3)
    ax.errorbar(mean_tau, TARGET, xerr=sd_tau, color="0.20", ecolor="0.20",
                marker="o", mfc=colour, mec="white", mew=1.2, ms=8,
                capsize=6, capthick=1.5, elinewidth=1.5, zorder=4)

    ax.text(.04, .95, "CL&P-scaled", transform=ax.transAxes, color=colour,
            fontsize=16, fontweight="bold", ha="left", va="top")
    ax.text(.04, .885, "310 K, five independent repeats", transform=ax.transAxes,
            fontsize=14, ha="left", va="top")
    ax.annotate(r"$C(t)=e^{-1}$", xy=(xmax, TARGET), xytext=(-8, 7),
                textcoords="offset points", ha="right", va="bottom", fontsize=14)
    ax.annotate("mean first crossing\n" + rf"$\tau_{{1/e}}={mean_tau:.1f}\pm{sd_tau:.1f}$ ps",
                xy=(mean_tau, TARGET), xytext=(128, .54), textcoords="data",
                arrowprops={"arrowstyle": "-", "lw": 1.2, "color": "0.25"},
                ha="left", va="bottom", fontsize=14)

    ax.set(xlim=(0, xmax), ylim=(0, 1.05), xlabel="Time (ps)",
           ylabel=r"Re-forming nearest-neighbour $C(t)$")
    ax.set_yticks(np.arange(0, 1.1, .2))
    ax.tick_params(width=1.5, length=7)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    fig.tight_layout()
    fig.savefig(OUTPUT, dpi=300, bbox_inches="tight", facecolor="white")


if __name__ == "__main__":
    main()
