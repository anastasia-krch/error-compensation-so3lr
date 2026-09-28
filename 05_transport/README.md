# Transport properties

`data/transport_vs_temperature.csv` contains the simulation values and uncertainties used for the self-diffusion, Green--Kubo conductivity, Nernst--Einstein conductivity, and Haven-ratio panels. Experimental reference data from Tokuda et al. are in `reference_data/`.

Transport coefficients and their uncertainties are estimated from velocity or charge-current time series with STACIE. `scripts/extract_transport_arrays.py` prepares these time series, `scripts/fit_stacie.py` performs the spectral fits, and `scripts/build_transport_table.py` writes the plotting table.

Run `notebooks/plot_transport.ipynb` to regenerate the PNG figures in `figures/`.
