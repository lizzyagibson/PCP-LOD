# Principle-Component-Pursuit
PCP for environmental health / epidemiology.

This version of PCP uses the non-convex rank approximation and the square root of the residual term in the objective function, along with distinct penalties for values less than the limit of detection (<LOD).

This repo includes:

0. `pcpr.R` file with PCP-LOD and CV functions
1. `functions.R` loads packages and includes some functions that I wrote
2. `NHANES` folder
    * `nhanes_cleaning.R` includes a function to clean the NHANES data and makes a figure
    * `nhanes_pops_50.Rmd` includes analysis and makes some figures
    * `Data` subfolder has original NHANES SAS file and R object with saved CV fits
3. `Sims` folder
    * `sim_grid.R` makes the simulated datasets
    * `sim_pca.R` runs PCA on all simulated datasets
    * `sim_ncvx_cv.R` cross-validated PCP-LOD on all simulated datasets
    * `sim_ncvx_pcp.R` runs PCP-LOD on all simulated datasets
    * `make_figures.R` makes figures for the manuscript
    * `Sim Data` subfolder has CV fits
    
    
    