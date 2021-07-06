# Principle-Component-Pursuit
PCP for environmental health / epidemiology.

This version of PCP uses the non-convex rank approximation and the square root of the residual term in the objective function, along with distinct penalties for values less than the limit of detection (<LOD).

You will need two packages, `pcpr` and `PCPhelpers` from the Columbia Prime GitHub page. To use these packages, clone the repos or download the .zip files. Locate the folders, unzip if applicable, and use the following code. Replace `path_to_folder` with your local path.

`install.packages("path_to_folder/pcpr", repos = NULL, type="source")`
`install.packages("path_to_folder/PCPhelpers", repos = NULL, type="source")`
`library(pcpr)`
`library(PCPhelpers)`

This repo includes:

1. `functions.R' loads packages and includes some functions that I wrote.
2. `NHANES` folder
    * `nhanes_cleaning.R` includes a function to clean the NHANES data and makes a figure
    * `nhanes_pops_50.Rmd` includes analysis and makes some figures
3. Sims
    * `sim_grid.R` makes the simulated datasets
    * `sim_pca.R` runs PCA on all simulated datasets
    * `sim_ncvx_cv.R` cross-validated PCP-LOD on all simulated datasets
    * `sim_ncvx_pcp.R` runs PCP-LOD on all simulated datasets
    * `make_figures.R` makes figures for the manuscript
    
    
    
    