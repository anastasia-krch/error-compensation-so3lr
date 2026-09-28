# Runbook

Run the shared extractor once per raw trajectory, then work through the paper
chapters in numerical order. Keep raw trajectories, reduced XYZ files, Lion
outputs, and binary transport arrays outside this repository.

1. `01_density`: extract density, build the compact table, execute the notebook.
2. `02_rdf`: make reduced nitrogen coordinates, run/pool RDFs, execute the notebook.
3. `03_ion_pair_correlations`: compute nearest-neighbour intermittent correlations and plot the compact summary.
4. `04_hydrogen_bonds`: calculate/pool H-bond RDFs and execute the plot notebook.
5. `05_transport`: extract charge-current and ion-velocity arrays externally, fit STACIE properties, export a compact table, and execute the transport notebook.
6. `06_dihedrals`: extract/pool the NTf2 dihedral histogram and execute the dihedral notebook.
7. `07_charges`: prepare the frozen clusters, run every method through the AMS driver (SO3LR through Engine ASE), compute the three charge caches, import the compact tables, and execute the charge notebook. The committed CSVs already reproduce the figures, so the notebook runs without AMS or a SO3LR installation.

The force-field argument is always explicit: `apple`, `so3lr`, `clp`, or
`clp-scaled`. The shared extraction code contains all topology maps.
