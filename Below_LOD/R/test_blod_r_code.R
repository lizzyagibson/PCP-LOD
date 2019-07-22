getwd()
source("./Below_LOD/R/function_blod.R")

library(R.matlab)
library(tidyverse)
options(scipen = 999)

# Fake matrix to try each function on

fake <- matrix(1:100, nrow = 10, byrow = TRUE)
fake
c <- 0.1

# prox_l1
R_l1 <- prox_l1(fake, c)

M_l1 <- readMat("./Below_LOD/MATLAB/LOD_demo_output/test_M_l1.mat") %>% 
  as.data.frame() %>% as_tibble() %>% as.matrix()

R_l1 - M_l1
# SAME!

#prox_nuclear
R_nuc <- prox_nuclear(fake, c)

M_nuc <- readMat("./Below_LOD/MATLAB/LOD_demo_output/test_M_nuc.mat") %>% 
  as.data.frame() %>% as_tibble() %>% as.matrix()

R_nuc[[1]] - M_nuc
# SAME!

# is_same
f1 <- matrix(1:100, nrow = 10, byrow = TRUE)
f2 <- matrix(1:100, nrow = 10, byrow = TRUE)
f3 <- matrix(5:104, nrow = 10, byrow = TRUE)
  
MAX_ITER <- 100
LOSS_THRESH <- 1e-5
SAME_THRESH <- 1e-4

is_same(SAME_THRESH, fake, f1, f3)
#SAME!

# loss_lod

# pcp_lod

pcp_lod(matrix(1:4, nrow = 2), 4/sqrt(10), 1, 0)

###################################################
###################################################

mixture <- readMat("./Data/mixtures_data.mat")

# Original data, none below LOD
mix <- as.data.frame(mixture) %>% as_tibble() %>% select(Al, As, Ba, bc, Br, Ca, Cl,
                                                                  Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
                                                                  Ti,  V, Zn)

mix <- as.matrix(mix[complete.cases(mix),])
D <- mix

m <- nrow(mix)
n <- ncol(mix)
lambda <- 4/sqrt(m)
mu <- 10
Delta <- 0

R_out <- pcp_lod(scale(mix), 4/sqrt(m), 10, 0)
R_S_lod0 <- R_out$S
R_L_lod0 <- R_out$L

# Read in MATLAB results
M_L_lod0 <- readMat("./Below_LOD/MATLAB/LOD_demo_output/lowrank_lod0.mat") %>% 
  as.data.frame() %>% as_tibble() %>% as.matrix()
M_S_lod0 <- readMat("./Below_LOD/MATLAB/LOD_demo_output/sparse_lod0.mat") %>% 
  as.data.frame() %>% as_tibble() %>% as.matrix()

# Compare with MATLAB results

## S matrices

colSums(R_S_lod0)
colSums(M_S_lod0)

## L matrices

colSums(R_L_lod0)
colSums(M_L_lod0)

# svd(L)
colSums(svd(R_L_lod0)$v)
colSums(svd(M_L_lod0)$v)

svd(R_L_lod0)$d
svd(M_L_lod0)$d





