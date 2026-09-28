# Environment-dependent SO3LR charges

This chapter examines SO3LR ionic charges, molecular dipoles, and ion-pair dissociation at fixed geometries.

`figures/charge_convergence.png` shows the mean ionic charge for neutral clusters around six spatially sampled central ion pairs. Shading is one standard error over those samples.

`figures/neutral_control.png` compares SO3LR and DFT dipoles using 50 ion-pair geometries and seven neutral molecules.

`figures/pes_scan.png` compares the frozen ion-pair dissociation curves from PBE0-D3(BJ)/TZ2P, APPLE&P, and SO3LR.

The plotted values are in `data/charge_convergence.csv`, `data/neutral_control.csv`, and `data/pes_scan.csv`. Run `notebooks/plot_charges.ipynb` to regenerate the three PNG figures.
