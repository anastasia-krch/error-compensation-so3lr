# [EMIM][NTf2] paper reproducibility repository

Start with [RUNBOOK.md](RUNBOOK.md). This repository contains analysis and
simulation setup code, supplied references, compact plot data, and executed
notebooks. Raw trajectories, reduced XYZ files, and binary transport arrays
stay outside it.

| Chapter | Content |
|---|---|
| `01_density` | Density versus temperature |
| `02_rdf` | Nitrogen RDFs |
| `03_ion_pair_correlations` | Intermittent nearest-neighbour correlations |
| `04_hydrogen_bonds` | H-bond RDFs |
| `05_transport` | Diffusion, conductivity, NE conductivity, and Haven ratio |
| `06_dihedrals` | NTf2 backbone dihedral distributions |
| `07_charges` | Partial charges: cluster convergence, dipole control, dissociation curve |
| `simulations` | AMS/APPLE&P and SO3LR NPT -> NVT -> NVE setup templates |

`analysis/scripts/extract_reduced_xyz.py` is the shared structural extractor.
All force-field maps are explicit in source code: APPLE&P uses its own raw
order; SO3LR, CL&P, and CL&P-scaled use the standard declared order.

Simulation templates are in `simulations/`. They are parameterised rather
than tied to a particular cluster or storage location.
