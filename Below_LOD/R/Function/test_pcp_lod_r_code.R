getwd()
source("./Below_LOD/R/function_pcp_lod.R")

library(R.matlab)
library(tidyverse)
options(scipen = 999)

# Fake matrix to try each function on
fake <- matrix(1:100, nrow = 10, byrow = TRUE)
c <- 0.1

# pcp_lod
pcp_lod(fake, 1/sqrt(10), 10, 0)
#SAME!

# prox_l1
prox_l1(fake, c)
# SAME!

#prox_nuclear
prox_nuclear(fake, c)
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
loss_lod(fake, f3, 0)
#SAME!

###################################################
###################################################

#real data
mixture <- readMat("./Data/mixtures_data.mat")

# Original data, none below LOD
mix <- as.data.frame(mixture) %>% as_tibble() %>% select(Al, As, Ba, bc, Br, Ca, Cl,
                                                                  Cr, Cu, Fe, K, Mn,  Ni,  Pb,  S,  Se,  Si,
                                                                  Ti,  V, Zn)

mix <- as.matrix(mix[complete.cases(mix),])
m <- nrow(mix)
n <- ncol(mix)

R_out <- pcp_lod(mix, 4/sqrt(m), 10, rep(0, 20))
R_S_lod0 <- R_out$S
R_L_lod0 <- R_out$L

# Read in MATLAB results
M_L_lod0 <- readMat("./Below_LOD/MATLAB/LOD_demo_output/lowrank_lod0.mat") %>% 
  as.data.frame() %>% as_tibble() %>% as.matrix()
M_S_lod0 <- readMat("./Below_LOD/MATLAB/LOD_demo_output/sparse_lod0.mat") %>% 
  as.data.frame() %>% as_tibble() %>% as.matrix()

# Compare with MATLAB results

## S matrices
head(R_S_lod0)
head(M_S_lod0)
R_S_lod0 - M_S_lod0

## L matrices
head(R_L_lod0)
head(M_L_lod0)

# svd(L)
colSums(svd(R_L_lod0)$v)
colSums(svd(M_L_lod0)$v)

svd(R_L_lod0)$d
svd(M_L_lod0)$d

### NMF
library(NMF)

results_nmf <- nmf(R_L_lod0, rank = 5)
sum(R_L_lod0 < 0)
R_L_lod0[R_L_lod0 < 0]

summary(results_nmf)
fit(results_nmf)
# W
basis(results_nmf)
# H
coef(results_nmf)
