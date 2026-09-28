# Hydrogen-bond radial distribution functions

This chapter compares the 300 K RDFs between all ring C--H donors of EMIM+ and O or F atoms of NTf2-.

`data/hbond_rdfs_300K.csv` contains the plotted simulation curves. The AIMD reference curves from Brela et al. are `reference_data/O-Hbond.csv` and `reference_data/O-Fbond.csv`. Run `notebooks/plot_hbond_rdfs_300K.ipynb` to recreate `figures/hbond_rdfs_300K.png`.

The analysis scripts accept force-field-specific reduced trajectories and produce the compact RDF table. Atom selections are defined in `analysis/scripts/extract_reduced_xyz.py`.
