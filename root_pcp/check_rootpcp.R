library(R.matlab)
library(tidyverse)
library(pcpr)
library(tictoc)
mixture <- readMat("/Users/lizzy/Principal.Component.Pursuit/Data/mixtures_data.mat")

mixture_data <- as.data.frame(mixture) %>% 
  as_tibble() %>% dplyr::select(Al, As, Ba, bc, Br, Ca, Cl,
                         Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
                         Ti,  V, Zn) %>% as.matrix()

X <- mixture_data[complete.cases(mixture_data),]
# X[1:5, 1:5]

m <- nrow(X)
n <- ncol(X)

tic()
root_out = pcp_lod(X, 1/sqrt(m), 10, 0)
toc()

tic()
root_out = root_pcp(X, 1/sqrt(m), 10)
toc()

tic()
root_nn = root_pcp_nonnegL(X, 1/sqrt(m), 10)
toc()

norm(root_out[[1]], "F")
norm(root_out[[2]], "F")

norm(root_nn[[1]], "F")
norm(root_nn[[2]], "F")
