# scripts/

* AnnualRecruitsSampledTable.R: Script from [Clownfish_persistence](https://github.com/pinskylab/Clownfish_persistence), used to estimate the number of recruits sampled at each site.

* BioPhysicalMatchingMatlabBode.txt: MATLAB code from Bode et al. 2019 PLOS Biology. Used as a guide for writing the likelihood function in R.
* Format-AddDest-Tables-For-Likelihood-Input.R: Script that creates destination site data (2022-01-20_AddDestTables.Rdata) used in pseudoR2.R, examine_LL_by_site.R bootstrap_data.ipynb
* Format-Input-Survey-Data.R: Code to format the genetic survey data to be input into the likelihood functions used in fit_data.ipynb
* Format-Simulation-Table-For-Likelihood-Input.R: Code to format the biophysical simulation data to be input into the likelihood functions used in fit_data.ipynb
* GenGausKernInt_sum0.5.R: Script to integrate kernels to 0.5 (used in cdf_solve.R for notebooks to get dispersal distances from kernels in notebooks fit_data.ipynb and bootstrap_data.ipynb
* LL_biophys.R: Full likelihood function used to fit dispersal kernels
* PropSampTable.R: Code to calculate the proportion of habitat sampled each year for input into the likelihood equation ProportionHabitatSampled.csv
* Seeded-Particle-Info.R: Script to calculate the total number of particles seeded by the model at each grid cell in each time point
* cdf_solve.R: Script to calculate the dispersal distances from kernel fits
* cdf_solve90.R: Script to calculate the distance at which 90% of dispersal occurs given kernel fit
* drifter_analysis.R: Plot drift tracks from 2017 deployment
* format_biophys_normalized_matrix.R: Script to make input parentage matrix for all biophysical results, normalizing based on the number of particles seeded at each site
* format_biophys_parentage_matrix.R: Script to make input parentage matrix for all biophysical results- CAI included with Other as "unknown"
* format_data.R: Format survey data and biophysical simulation data similarly to input into the kernel fitting likelihood function
* format_genetic_kernel_parentage_matrix.R: Format the genetic parentage data for likelihood equation
* format_genetic_parentage_matrix.R: Unclear, looks to be the same as the above file. I think fine to delete
* ll_kt_both_bbmle.R: Function to estimate maximum likelihood k and theta using bbmle package
* ll_kt_both_grid_search.R: Function to estimate maximum likelihood k and theta using brute force grid search, just to validate bbmle results in testing phase of analysis
* plotEKE_FigS2.R: Make a plot of eddy kinetic energy by site by year and monsoon. Becomes Fig. S2 in the manuscript
* pseudoR2.R: calculate McFadden's pseudo-R2 from the log likelihoods of the parentage and connectivity simulations. A subset of the code in notebooks/bootstrap_data.ipynb
* roms_vertices.R: Code to join data in ROMS vertices input files and QGIS output files to create master tables mapping surveyed sites to ROMs model locations
