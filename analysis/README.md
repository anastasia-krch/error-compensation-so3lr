# Shared structural extraction

`scripts/extract_reduced_xyz.py` reads each raw AMS RKF or SO3LR HDF5
trajectory once and writes external reduced XYZ files for the nitrogen RDF,
ion-pair, hydrogen-bond, and dihedral chapters. Raw trajectories and reduced
XYZ files are not publication-repository data.

Choose the topology explicitly with `--force-field apple`, `so3lr`, `clp`, or
`clp-scaled`. The atom maps are defined in the extractor itself; APPLE&P has
its own raw ordering and the remaining force fields share the standard map.

```bash
amspython analysis/scripts/extract_reduced_xyz.py --source rkf --force-field apple \
  --rkf /path/to/ams.rkf --out-prefix /external/apple_300K_seed1
```

```bash
python analysis/scripts/extract_reduced_xyz.py --source hdf5 --force-field so3lr \
  --hdf5 /path/to/nve_part*.hdf5 --out-prefix /external/so3lr_300K_seed0
```
