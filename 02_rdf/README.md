# Nitrogen RDFs

This chapter plots the 300 K anion--anion, cation-midpoint--cation-midpoint,
and cation-midpoint--anion nitrogen RDFs. The shared extractor makes external
`*_nmid_anionN.xyz` files with central NTf2 N followed by EMIM ring-N
midpoints, using the explicit force-field map selected at extraction time.

Run `scripts/run_lionanalysis_all_nitrogen_rdfs.py`, then
`scripts/pool_all_nitrogen_rdfs.py`, and execute
`notebooks/plot_nitrogen_rdfs_300K.ipynb`. The committed CSV and figure are
the compact plotting artefacts; trajectories and Lion output remain external.
