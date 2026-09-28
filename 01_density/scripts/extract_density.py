#!/usr/bin/env python3
"""Extract NPT density from AMS (APPLE&P) or SO3LR MD log files.

Neither log stores density directly; it is computed from the simulation-cell
volume in each frame and the (constant) total system mass:

    density [g/cm^3] = mass_g / (volume_A3 * 1e-24)

Two log formats are supported:

  * AMS (`ams.log`): whitespace columns after a `<date> <time>` prefix, i.e.
        Step  Time(fs)  Temp(K)  E_pot(au)  Pressure(MPa)  Volume(A^3)

  * SO3LR (`npt_jax.log`): `... SO3LR - INFO -` lines with columns
        Step  E  KE  PE  H  Temp(K)  (lx,ly,lz)  time/step

The averaging window is the last fraction (default 0.5 = "last half") of the
trajectory by time, and the reported error bar is the standard error of the mean
from block averaging (default 10 blocks), which accounts for autocorrelation of
NPT frames.

Depends only on numpy so it is easy to copy to another machine.

Examples
--------
    python extract_density.py --engine so3lr --log 290K/npt_jax.log
    python extract_density.py --engine ams --log npt_fromkrf_290K.results/ams.log \
        --mass-amu 97827.765
    python extract_density.py --engine ams --log .../ams.log --xyz .../start.xyz
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
# Physical constants
# --------------------------------------------------------------------------- #
N_AVOGADRO = 6.02214076e23
AMU_TO_G = 1.0 / N_AVOGADRO          # 1 amu in grams (== 1.66053906660e-24 g)
A3_TO_CM3 = 1.0e-24

# Atomic masses (g/mol) for the elements present in EMIM / NTf2.
ATOMIC_MASS = {
    "H": 1.00794,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984032,
    "S": 32.065,
}


# --------------------------------------------------------------------------- #
# Mass helpers
# --------------------------------------------------------------------------- #
def total_mass_amu_from_xyz(xyz_path: Path) -> float:
    """Sum atomic masses (amu) of all atoms in an XYZ file."""
    with Path(xyz_path).open() as handle:
        n_atoms = int(handle.readline().strip())
        handle.readline()  # comment line
        species = [handle.readline().split()[0] for _ in range(n_atoms)]

    missing = sorted(set(species) - set(ATOMIC_MASS))
    if missing:
        raise ValueError(f"Missing atomic masses for elements: {missing}")
    return sum(ATOMIC_MASS[element] for element in species)


# --------------------------------------------------------------------------- #
# Log parsers -> structured arrays
#   Each returns a dict of numpy arrays with keys:
#     step, time_ps, temp_K, vol_A3, density_g_cm3
#   AMS additionally provides press_MPa.
# --------------------------------------------------------------------------- #
def parse_ams_log(path: Path, mass_amu: float) -> dict:
    """Parse an AMS `ams.log` MD run.

    Data rows look like:
        <May12-2026> <16:47:11>  0  0.00  289.  -37.88667  -2.285  106026.2
    columns after the two `<...>` timestamp tokens:
        Step  Time(fs)  Temp(K)  E_pot(au)  Pressure(MPa)  Volume(A^3)
    Header and non-data lines are skipped.
    """
    mass_g = mass_amu * AMU_TO_G
    steps, time_fs, temp, press, vol = [], [], [], [], []

    with Path(path).open() as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 8:
                continue
            # Strip the two `<date> <time>` tokens when present.
            data = parts[2:8] if parts[0].startswith("<") and parts[1].startswith("<") else parts[:6]
            if len(data) != 6 or not data[0].isdigit():
                continue
            try:
                s = int(data[0])
                t_fs, t_k, _e_au, p_mpa, v_a3 = map(float, data[1:])
            except ValueError:
                continue
            steps.append(s)
            time_fs.append(t_fs)
            temp.append(t_k)
            press.append(p_mpa)
            vol.append(v_a3)

    if not steps:
        raise ValueError(f"No MD data rows found in AMS log: {path}")

    vol = np.asarray(vol, dtype=float)
    return {
        "step": np.asarray(steps, dtype=np.int64),
        "time_ps": np.asarray(time_fs, dtype=float) / 1000.0,
        "temp_K": np.asarray(temp, dtype=float),
        "press_MPa": np.asarray(press, dtype=float),
        "vol_A3": vol,
        "density_g_cm3": mass_g / (vol * A3_TO_CM3),
    }


_SO3LR_TEMP_RE = re.compile(r"Temperature:\s+([0-9.]+)\s*K")


def parse_so3lr_log(path: Path, mass_amu: float, dt_fs: float = 0.5) -> dict:
    """Parse a SO3LR `npt_jax.log` NPT run.

    Data lines look like (after the logging prefix):
        ... SO3LR - INFO - 100  -12182.99  319.79  -12502.79  -12184.54  291.1  (46.037,46.895,45.005)  1.82e-01
    columns: Step  E  KE  PE  H  Temp(K)  (lx,ly,lz)  time/step
    """
    mass_g = mass_amu * AMU_TO_G
    steps, temp, vol = [], [], []

    with Path(path).open() as handle:
        for line in handle:
            if "SO3LR - INFO -" not in line:
                continue
            payload = line.split("SO3LR - INFO -", 1)[1].strip()
            fields = payload.split()
            if len(fields) < 8 or not fields[0].isdigit():
                continue
            box_field = fields[6]
            if not (box_field.startswith("(") and box_field.endswith(")")):
                continue
            try:
                lx, ly, lz = (float(v) for v in box_field.strip("()").split(","))
                s = int(fields[0])
                t_k = float(fields[5])
            except ValueError:
                continue
            steps.append(s)
            temp.append(t_k)
            vol.append(lx * ly * lz)

    if not steps:
        raise ValueError(f"No MD data rows found in SO3LR log: {path}")

    steps_arr = np.asarray(steps, dtype=np.int64)
    vol = np.asarray(vol, dtype=float)
    return {
        "step": steps_arr,
        "time_ps": steps_arr.astype(float) * dt_fs / 1000.0,
        "temp_K": np.asarray(temp, dtype=float),
        "press_MPa": None,  # not reported in SO3LR log
        "vol_A3": vol,
        "density_g_cm3": mass_g / (vol * A3_TO_CM3),
    }


def read_so3lr_target_temp(path: Path):
    """Return the target temperature (K) from a SO3LR log header, or None."""
    with Path(path).open() as handle:
        for line in handle:
            m = _SO3LR_TEMP_RE.search(line)
            if m:
                return float(m.group(1))
            if "SO3LR - INFO - 0\t" in line or "SO3LR - INFO - 0 " in line:
                break  # past the header
    return None


# --------------------------------------------------------------------------- #
# Statistics helpers
# --------------------------------------------------------------------------- #
def last_fraction_mask(time_ps: np.ndarray, fraction: float) -> np.ndarray:
    """Boolean mask selecting the last `fraction` of the trajectory by time."""
    if not 0.0 < fraction <= 1.0:
        raise ValueError("--fraction must be in (0, 1]")
    t0, t1 = float(time_ps[0]), float(time_ps[-1])
    cutoff = t1 - fraction * (t1 - t0)
    return time_ps >= cutoff


def block_sem(values: np.ndarray, n_blocks: int) -> float:
    """Standard error of the mean from `n_blocks` equal-size blocks.

    SEM = std(block_means, ddof=1) / sqrt(n_blocks).
    Trailing frames that do not fill a whole block are dropped.
    """
    values = np.asarray(values, dtype=float)
    block_size = len(values) // n_blocks
    if block_size < 1:
        raise ValueError(
            f"Not enough frames ({len(values)}) for {n_blocks} blocks"
        )
    trimmed = values[: block_size * n_blocks].reshape(n_blocks, block_size)
    block_means = trimmed.mean(axis=1)
    return float(np.std(block_means, ddof=1) / math.sqrt(n_blocks))


# --------------------------------------------------------------------------- #
# High-level API
# --------------------------------------------------------------------------- #
def analyze(engine: str, log: Path, mass_amu: float,
            fraction: float = 0.5, n_blocks: int = 10) -> dict:
    """Parse a log and return a summary dict for the last-`fraction` window."""
    if engine == "ams":
        data = parse_ams_log(log, mass_amu)
    elif engine == "so3lr":
        data = parse_so3lr_log(log, mass_amu)
    else:
        raise ValueError(f"Unknown engine: {engine!r}")

    mask = last_fraction_mask(data["time_ps"], fraction)
    rho = data["density_g_cm3"][mask]
    temp = data["temp_K"][mask]
    win_time = data["time_ps"][mask]

    press_mean = ""
    if data["press_MPa"] is not None:
        press_mean = float(np.mean(data["press_MPa"][mask]))

    return {
        "engine": engine,
        "mass_amu": mass_amu,
        "n_frames_total": int(len(data["time_ps"])),
        "n_frames_window": int(rho.size),
        "window_start_ps": float(win_time[0]),
        "window_end_ps": float(win_time[-1]),
        "density_mean_g_cm3": float(np.mean(rho)),
        "density_sem_g_cm3": block_sem(rho, n_blocks),
        "density_std_g_cm3": float(np.std(rho, ddof=1)),
        "mean_temp_K": float(np.mean(temp)),
        "mean_press_MPa": press_mean,
        "_data": data,          # full arrays, for optional CSV export
        "_mask": mask,
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _infer_ams_temp(log: Path):
    m = re.search(r"npt_fromkrf_(\d+)K", str(log))
    return float(m.group(1)) if m else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Extract NPT density (last-half mean +/- block SEM) from an "
                    "AMS or SO3LR MD log.")
    ap.add_argument("--engine", required=True, choices=["ams", "so3lr"])
    ap.add_argument("--log", required=True, type=Path, help="path to the log file")
    ap.add_argument("--mass-amu", type=float, default=None,
                    help="total system mass in amu (overrides --xyz)")
    ap.add_argument("--xyz", type=Path, default=None,
                    help="XYZ file to derive total mass from")
    ap.add_argument("--fraction", type=float, default=0.5,
                    help="last fraction of trajectory to average (default 0.5)")
    ap.add_argument("--blocks", type=int, default=10,
                    help="number of blocks for SEM (default 10)")
    ap.add_argument("--target-temp", type=float, default=None,
                    help="label temperature (K); auto-detected if omitted")
    ap.add_argument("--out-csv", type=Path, default=None,
                    help="optional per-frame timeseries CSV to write")
    args = ap.parse_args(argv)

    # Resolve total mass.
    if args.mass_amu is not None:
        mass_amu = args.mass_amu
    else:
        xyz = args.xyz
        if xyz is None and args.engine == "so3lr":
            xyz = args.log.parent / "start.xyz"  # sibling default for SO3LR
        if xyz is None or not Path(xyz).exists():
            ap.error("provide --mass-amu, or --xyz pointing to a start.xyz")
        mass_amu = total_mass_amu_from_xyz(xyz)

    res = analyze(args.engine, args.log, mass_amu,
                  fraction=args.fraction, n_blocks=args.blocks)

    # Resolve a temperature label for the printout.
    target = args.target_temp
    if target is None:
        target = (read_so3lr_target_temp(args.log) if args.engine == "so3lr"
                  else _infer_ams_temp(args.log))
    target_str = f"{target:.0f}K" if target is not None else "?"

    print(f"engine={res['engine']}  T_target={target_str}  "
          f"mass={mass_amu:.3f} amu")
    print(f"frames: {res['n_frames_window']} of {res['n_frames_total']} "
          f"in window [{res['window_start_ps']:.1f}, "
          f"{res['window_end_ps']:.1f}] ps (last {args.fraction:g})")
    print(f"density = {res['density_mean_g_cm3']:.5f} "
          f"+/- {res['density_sem_g_cm3']:.5f} g/cm^3  "
          f"(block SEM, {args.blocks} blocks; "
          f"frame std {res['density_std_g_cm3']:.5f})")
    print(f"mean T in window = {res['mean_temp_K']:.2f} K")

    if args.out_csv is not None:
        data = res["_data"]
        cols = [data["step"], data["time_ps"], data["temp_K"],
                data["vol_A3"], data["density_g_cm3"]]
        header = "Step,Time_ps,Temp_K,Vol_A3,Density_g_cm3"
        fmt = ["%d", "%.6f", "%.6f", "%.6f", "%.8f"]
        np.savetxt(args.out_csv, np.column_stack(cols),
                   delimiter=",", header=header, comments="", fmt=fmt)
        print(f"wrote per-frame timeseries: {args.out_csv}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
