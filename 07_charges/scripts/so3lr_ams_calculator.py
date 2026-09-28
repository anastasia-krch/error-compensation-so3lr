#!/usr/bin/env python3
"""Charge-aware SO3LR calculator for AMS's ``Engine ASE``.

The stock SO3LR ASE calculator exposes energy and forces only.  This adapter asks the
network for its intermediate ``partial_charges`` output and forwards it, and the
charge-derived dipole, through the ASE/AMS pipe.  AMS therefore owns the calculation
and its RKF/log artefacts just as it does for ADF and APPLE&P.

The file is deliberately self-contained: AMS imports it in the isolated SO3LR uv
environment, where sibling modules from this repository are not on ``sys.path``.
"""
from __future__ import annotations

import numpy as np

from ase.calculators.calculator import Calculator, all_changes
from scm.amspipe import AMSExternalCapabilities


class So3lrChargeCalculator(Calculator):
    """SO3LR with energy, forces, atomic charges, and a molecular dipole."""

    implemented_properties = ["energy", "forces", "stress", "free_energy", "hessian",
                              "charges", "dipole"]

    def __init__(self, *, lr_cutoff: float = 12.0,
                 dispersion_energy_cutoff_lr_damping: float = 2.0,
                 calculate_stress: bool = False, dtype: str = "float64", **kwargs) -> None:
        super().__init__()
        from so3lr import So3lrCalculator

        precision = {"float32": np.float32, "float64": np.float64}.get(str(dtype))
        if precision is None:
            raise ValueError("dtype must be 'float32' or 'float64'")
        # ``has_aux`` keeps the network's intermediate quantities in the underlying
        # calculator result.  ``partial_charges`` is exactly the quantity used by the
        # original direct-JAX analysis, not an AMS charge partitioning scheme.
        self._calculator = So3lrCalculator(
            lr_cutoff=lr_cutoff,
            dispersion_energy_cutoff_lr_damping=dispersion_energy_cutoff_lr_damping,
            calculate_stress=calculate_stress,
            dtype=precision,
            observables=["partial_charges"],
            has_aux=True,
            **kwargs)
        capabilities = AMSExternalCapabilities(
            charges=True,
            periodicity=True,
            charged_periodic_cell=True)
        capabilities.apply_implemented_properties(self.implemented_properties)
        self.ams_capabilities = capabilities

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self._calculator.calculate(atoms, properties, system_changes)
        result = self._calculator.results
        charges = np.asarray(result["aux"]["partial_charges"], dtype=float).reshape(-1)
        # ASE stores dipoles in e Angstrom.  The AMS pipe converts that to e bohr.
        dipole = np.einsum("i,ij->j", charges, atoms.get_positions())
        self.results = {
            key: value for key, value in result.items() if key != "aux"
        }
        self.results["charges"] = charges
        self.results["dipole"] = dipole


def get_calculator(**kwargs):
    """AMS-required factory function."""
    return So3lrChargeCalculator(**kwargs)
