#!/usr/bin/env python3
"""Fit diffusion and Green--Kubo conductivity from one external array set."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.constants import Boltzmann, elementary_charge


def fit(sequences, timestep_s, prefactor, symbol, unit):
    from stacie import ExpPolyModel, UnitConfig, compute_spectrum, estimate_acint
    spectrum = compute_spectrum(sequences, prefactors=prefactor, timestep=timestep_s, include_zero_freq=False)
    units = UnitConfig(acint_symbol=symbol, acint_unit_str=unit, acint_fmt='.3e',
                       time_unit=1e-12, time_unit_str='ps', freq_unit=1e12, freq_unit_str='THz')
    result = estimate_acint(spectrum, ExpPolyModel([0, 1, 2]), verbose=True, uc=units)
    return {"value": float(result.acint), "std": float(result.acint_std),
            "fcut_Hz": float(result.fcut), "ncut": int(result.ncut), "neff": float(result.neff),
            "cost_zscore": float(result.props['cost_zscore']), "criterion_zscore": float(result.props['criterion_zscore'])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--temperature-K", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    meta = json.loads(args.metadata.read_text())
    velocities = np.load(meta["ionvels"], mmap_mode='r')
    current = np.load(meta["chargecurrent"], mmap_mode='r')
    timestep_s = float(meta["saved_dt_ps"]) * 1e-12
    factor = float(meta["velocity_to_m_s"])
    cation = np.ascontiguousarray(np.moveaxis(np.asarray(velocities[:, :250]) * factor, 0, -1).reshape(-1, len(velocities)))
    anion = np.ascontiguousarray(np.moveaxis(np.asarray(velocities[:, 250:]) * factor, 0, -1).reshape(-1, len(velocities)))
    current_si = np.ascontiguousarray((np.asarray(current) * factor * elementary_charge).T)
    volume_m3 = float(meta["volume_A3"]) * 1e-30
    result = {"temperature_K": args.temperature_K, "volume_A3": meta["volume_A3"],
              "cation_charge_e": meta["q_cation_e"], "anion_charge_e": meta["q_anion_e"],
              "diffusion_cation_m2_s": fit([cation], timestep_s, 1.0, 'D', 'm$^2$/s'),
              "diffusion_anion_m2_s": fit([anion], timestep_s, 1.0, 'D', 'm$^2$/s'),
              "conductivity_S_m": fit([current_si], timestep_s, 1 / (volume_m3 * Boltzmann * args.temperature_K), r'\\sigma', 'S/m')}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(f"wrote {args.output}")


if __name__ == '__main__': main()
