# Camotes Sea ROMS model - biophysical simulations of larval Amphiprion clarkii connectivity
Katrina Catalano's project calculating predicted dispersal using empirical oceanographic measurements and ROMS output for the Camotes Sea, Philippines, then comparing this to genetic estimates of connectivity.

Catalano, K., E. Drenkard, E. Curchitser, A. Dedrick, M. Stuart, H. Montes Jr., and M. L. Pinsky (_in press_). The contribution of nearshore oceanography to temporal variation in larval dispersal. **Ecology**

[![DOI](https://zenodo.org/badge/138894621.svg)](https://zenodo.org/doi/10.5281/zenodo.12744525)

This repository contains:
- notebooks/ Jupyter notebooks with R code used for analysis, in order of workflow
  - camotes_vertices_water.ipynb: Assign known source and destination sites in the study region to horizontal grid cells of the Camotes-ROMS model.  
  - process_ROMS.ipynb: Read in the Camotes-ROMS connectivity simulation output, add available physical metadata for sites, and make summary plots of simulation results. 
  - bootstrap_data.ipynb: Read in the Camotes-ROMS connectivity simulation output with metadata, format to fit dispersal kernels and make source normalized biophysical connectivity matrices for event matching. Also contains code for bootstrapping the simulation dispersal kernel fits.
  - fit_data.ipynb: Maximum likelihood estimation of ROMS dispersal kernels across time scales.
  - process_ROMS_kernels.ipynb: Read in the MLE dispersal kernel fits, add confidence intervals, calculate summary statistics, make final plots, and look at some potential climate correlations.
- script_output/ Content generated from code in notebooks/
- empirical_data/ Data from field (genetic) observations of dispersal, used for comparison to Camotes-ROMS connectivity simulations.
- ROMS/ Data input and output from ROMS model.
- scripts/ Raw R code versions of content in Jupyter notebooks, scripts called by notebooks, and a script to make Fig. S2 (eddy kinetic energy).

Code was run in R with Jupyter notebooks on a laptop and on a Centos Linux workstation.

Please contact mpinsky[at]ucsc.edu with questions.

## Data Accessibility
For data not on this git repo:

Pinsky, M., Stuart, M. (2023) Temperature loggers (HOBO) placed in two locations off the coast of the West coast of Leyte, the Philippines , 2012-2019. Biological and Chemical Oceanography Data Management Office (BCO-DMO). (Version 2) Version Date 2023-03-06. https://doi.org/10.26008/1912/bco-dmo.862415.2 

Pinsky, M. L.; Drenkard, E. J.; Curchitser, E. (2023): CamotesROMS surface currents. figshare. Dataset. https://doi.org/10.6084/m9.figshare.22827059.v1
