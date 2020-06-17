## 
# Check LOD dimensions
# Matrix vs. Vector vs. Scalar
##

library(tidyverse)
library(R.matlab)
source("./Below_LOD/R/function_pcp_lod.R")
options(scipen = 999)

# Read in Boston air pollution data
mixture <- readMat("./Data/mixtures_data.mat")

mixture_data <- as.data.frame(mixture) %>% as_tibble() %>% 
  select(Al, As, Ba, bc, Br, Ca, Cl,
         Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
         Ti,  V, Zn) %>% 
  drop_na()

mixture_data

# Create version with 10% lowest values for each variable as below the LOD
mix_data_lod_10 <- mixture_data %>% 
  mutate_all(~ifelse(. <= quantile(., probs = .10), -1, .)) %>% as.matrix()

# LOD quantiles

# vector
vector10 <- mixture_data %>% 
  summarise_all(quantile, probs = .10) %>% as_vector()

# matrix
tf_lod_10 <- mix_data_lod_10 %>% 
  as_tibble() %>% 
  mutate_all(~ifelse(. == -1, TRUE, FALSE)) %>% 
  as.matrix()
# T/F dataset, T = <LOD
matrix10  <- as.matrix(tf_lod_10*mixture_data)

if (is.vector(vector10)) {
  tf = ifelse(mix_data_lod_10 < 0, TRUE, FALSE)
  LOD = t(t(tf) * vector10)
}

all.equal(tf, tf_lod_10)
all.equal(matrix10, LOD)

norm(matrix10 - LOD, 'F')

apply(matrix10, 2, max) - vector10
apply(LOD, 2, max) - vector10

# to make vector LOD a matrix
apply(t(t(tf_lod_10) * vector10), 2, max)
  
# save to use in matlab
#write_csv(as_tibble(matrix10), "./Below_LOD/R/BLOD_airpol_data/mix_data_lod_10_matrixlod.csv")

# Run PCPLOD on 3 LODs
lambda <- 1/sqrt(nrow(mix_data_lod_10))

# scalar LOD
out_s <- pcp_lod(mix_data_lod_10, lambda, 10, 0.01)
r_low_s <- out_s[[1]]
r_sparse_s <- out_s[[2]]

# vector LOD
out_v <- pcp_lod(mix_data_lod_10, lambda, 10, vector10)
r_low_v <- out_v[[1]]
r_sparse_v <- out_v[[2]]

# matrix LOD
out_m <- pcp_lod(mix_data_lod_10, lambda, 10, matrix10)
r_low_m <- out_m[[1]]
r_sparse_m <- out_m[[2]]

# Load matlab results

mat_low_s <- readMat("./Below_LOD/MATLAB/LOD_demo_output/lowrank_lod_10s.mat")[[1]]
mat_low_v <- readMat("./Below_LOD/MATLAB/LOD_demo_output/lowrank_lod10v.mat")[[1]]
mat_low_m <- readMat("./Below_LOD/MATLAB/LOD_demo_output/lowrank_lod_10m.mat")[[1]]

mat_sparse_s <- readMat("./Below_LOD/MATLAB/LOD_demo_output/sparse_lod_10s.mat")[[1]]
mat_sparse_v <- readMat("./Below_LOD/MATLAB/LOD_demo_output/sparse_lod10v.mat")[[1]]
mat_sparse_m <- readMat("./Below_LOD/MATLAB/LOD_demo_output/sparse_lod_10m.mat")[[1]]

# Compare

# scalar LOD
norm(r_low_s    - mat_low_s,    type = "F")
norm(r_sparse_s - mat_sparse_s, type = "F")

# vector LOD
norm(r_low_v    - mat_low_v,    type = "F")
norm(r_sparse_v - mat_sparse_v, type = "F")

# matrix LOD
norm(r_low_m    - mat_low_m,    type = "F")
norm(r_sparse_m - mat_sparse_m, type = "F")
