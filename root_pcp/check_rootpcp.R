library(R.matlab)
library(tidyverse)
library(pcpr)
library(tictoc)

mixture <- readMat("/Users/lizzy/OneDrive - cumc.columbia.edu/Principal.Component.Pursuit/Data/mixtures_data.mat")

mixture_data <- as.data.frame(mixture) %>% 
  as_tibble() %>% dplyr::select(Al, As, Ba, bc, Br, Ca, Cl,
                         Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
                         Ti,  V, Zn) %>% as.matrix()

X <- mixture_data[complete.cases(mixture_data),]
# X[1:5, 1:5]

m <- nrow(X)
n <- ncol(X)
lambda <- 1/sqrt(m)
mu <- 1

## With missing values!
Xmissing = X
Xmissing[1:1000,1:5] = NA

## Run all model versions
#tic()
# lod_out = pcp_lod(X, lambda, mu, 0)
# #"Converged in 69 iterations."
# #toc()
# norm(lod_out[[1]], "F")
# norm(lod_out[[2]], "F")

# tic()
# root_out = root_pcp_noncvx_na(Xmissing, lambda, mu, 5, verbose = TRUE)
# toc()
# # converged in 145
# norm(root_out[[1]], "F")
# norm(root_out[[2]], "F")

tic()
root_out = root_pcp_noncvx_nonnegL_na(Xmissing, lambda, mu, 5, verbose = TRUE)
toc()
# converged in 287
# 3.7 sec elapsed
norm(root_out[[1]], "F")
norm(root_out[[2]], "F")

# tic()
# root_out = root_pcp_noncvx(X, lambda, mu, 5)
# #"Converged in 837 iterations."
# toc()
# # 14.049 sec elapsed
# norm(root_out[[1]], "F")
# norm(root_out[[2]], "F")

# #tic()
# root_nn = root_pcp_nonnegL(X, lambda, mu)
# #"Converged in 646 iterations."
# #toc()
# norm(root_nn[[1]], "F")
# norm(root_nn[[2]], "F")
# 
# tic()
# nan_out = root_pcp_with_nan(X, lambda, mu)
# #"Converged in 837 iterations."
# toc()
# # 16.649 sec elapsed
# norm(nan_out[[1]], "F")
# norm(nan_out[[2]], "F")
# 
# tic()
# nan_out = root_pcp_na(Xmissing, lambda, mu)
# #"Converged in 915 iterations."
# toc()
# # 18.347 sec elapsed
# norm(nan_out[[1]], "F")
# norm(nan_out[[2]], "F")
# 
# tic()
# nan_nn_out = root_pcp_na_nonnegL(Xmissing, lambda, mu)
# #"Converged in 1016 iterations."
# toc()
# # 18.347 sec elapsed
# norm(nan_nn_out[[1]], "F")
# norm(nan_nn_out[[2]], "F")
# 
# tic()
# nan_nn_lod_out = root_pcp_na_nonnegL_LOD(Xmissing, lambda, mu, 0)
# #"Converged in 632 iterations."
# toc()
# # 18.347 sec elapsed
# norm(nan_nn_lod_out[[1]], "F")
# norm(nan_nn_lod_out[[2]], "F")
