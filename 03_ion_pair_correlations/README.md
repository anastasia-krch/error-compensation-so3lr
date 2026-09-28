# Ion-pair residence lifetimes

Ion pairs are defined by cation-centred nearest-neighbour identity, without a distance cutoff. Pair re-formation is allowed. The reported lifetime is the first crossing of the intermittent correlation function at \(C(t)=e^{-1}\); see `METHODS.md`.

`data/nearest_neighbor_tau_1e_all_forcefields.csv` contains the plotted means and seed standard deviations. The corresponding per-seed values, trajectory windows, and endpoint correlations are in `data/nearest_neighbor_tau_1e_all_forcefields_by_seed.csv`.

The manuscript figures are `figures/nearest_neighbor_lifetime_first_crossing_example.png` and `figures/nearest_neighbor_ionpair_tau_1e_all_forcefields_logy.png`. Use `scripts/compute_nearest_neighbour.py` for individual reduced trajectories and `notebooks/plot_ion_pair_lifetimes.ipynb` to regenerate the temperature comparison.
