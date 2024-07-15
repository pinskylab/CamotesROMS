notebooks/

Working notebooks for the main analysis, ordered according to analysis flow:

* process_ROMS.ipynb: Convert the NetCDF ROMS output files into tidy data tables with dates formatted for manipulation in R
* fit_data.ipynb: Fit biophysical simulation dispersal kernels to the ROMS data
* bootstrap_data.ipynb: Bootstrap the simulation data and build ensemble simulation kernel fits to compare with genetic data kernel fits
* process_ROMS_kernels.ipynb: Summarize the biophysical kernel fits and compare to the genetic kernel fits. Create plots for visualizing results.

Other notebooks used in exploratory phase of analysis:
* drifter_tracks.ipynb: Visualize drifter tracks from field deployment in 2017
* velocity_fields.ipynb: Processes surface velocity fields from the ROMS model and visualizes. Based on out-of-date data as of May 2023.
