# Density

This chapter contains the density-versus-temperature comparison for systems of 250 [EMIM][NTf2] ion pairs.

`data/density_vs_temperature.csv` contains the plotted simulation values and `reference_data/tokuda_density.csv` contains the experimental reference. Each simulated density is the mean over the final half of the NPT trajectory; its uncertainty is the standard error of ten consecutive block means.

`scripts/extract_density.py` calculates one value from an AMS or SO3LR NPT log. `scripts/build_density_table.py` applies the same calculation to a manifest of simulations and writes the plotting table.

Run `notebooks/plot_density.ipynb` to recreate `figures/density.png`.
