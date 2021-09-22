# Camotes Sea ROMS model - biophysical simulations of larval Amphiprion clarkii connectivity
Katrina Catalano's project calculating predicted dispersal using empirical oceanographic measurements and ROMS output for the Camotes Sea, Philippines, then comparing this to genetic estimates of connectivity.

This repository contains:
- notebooks/ Jupyter notebooks with R code used for analysis, in order of workflow
  - camotes_vertices_water.ipynb: Assign known source and destination sites in the study region to horizontal grid cells of the Camotes-ROMS model.  
  - process_ROMS.ipynb: Read in the Camotes-ROMS connectivity simulation output, add available physical metadata for sites, and make summary plots of simulation results. 
  - format_data.ipynb: Read in the Camotes-ROMS connectivity simulation output with metadata, format to fit dispersal kernels and make source normalized biophysical connectivity matrices for event matching. Also contains code for bootstrapping the simulation dispersal kernel fits.
  - fit_data.ipynb: Maximum likelihood estimation of ROMS dispersal kernels across time scales.
  - process_ROMS_kernels.ipynb: Read in the MLE dispersal kernel fits, add confidence intervals, calculate summary statistics, make final plots, and look at some potential climate correlations.
- script_output/ Content generated from code in notebooks/
- empirical_data/ Data from field (genetic) observations of dispersal, used for comparison to Camotes-ROMS connectivity simulations.
- ROMS/ Data input and output from ROMS model.
- scripts/ Raw R code verisions of content in notebooks/

Please contact kat.catalano[at]rutgers.edu with questions.
