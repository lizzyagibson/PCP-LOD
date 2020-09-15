#!/bin/bash
#$ -cwd -S /bin/bash
#$ -l mem=4G
#$ -l time=:1000:
#$ -M eag2186@cumc.columbia.edu

$MODULESHOME/init/bash
module load R/3.6.0

clear

R CMD BATCH --no-save parameter_grid_boston.R

