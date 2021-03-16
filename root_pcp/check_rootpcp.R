library(R.matlab)
library(tidyverse)
library(pcpr)
library(tictoc)

mixture <- readMat("./Data/mixtures_data.mat")

mixture_data <- as.data.frame(mixture) %>% 
  as_tibble() %>% dplyr::select(Al, As, Ba, bc, Br, Ca, Cl,
                         Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
                         Ti,  V, Zn) %>% as.matrix()

X <- mixture_data[complete.cases(mixture_data),]

m <- nrow(X)
p <- ncol(X)
lambda <- 1/sqrt(m)
mu <- sqrt(p/2)

## With missing values!
Xmissing = X
Xmissing[1:1000,1:5] = NA

## Run all model versions

# 9 rootPCP functions, so far

# 9
# Does not converge...
root_out = root_pcp_noncvx_nonnegL_na(Xmissing, lambda, mu, 5, verbose = TRUE)
norm(root_out[[1]], "F")
norm(root_out[[2]], "F")

# 8
# SAME
root_out2 = root_pcp_noncvx_na(Xmissing, lambda, mu, 5, verbose = TRUE)
# converged in 4103
norm(root_out2[[1]], "F")
norm(root_out2[[2]], "F")

# 7
# SAME
root_nn = root_pcp_noncvx_nonneg(X, lambda, mu, 5, verbose = TRUE)
#"Converged in 235 iterations."
norm(root_nn[[1]], "F")
norm(root_nn[[2]], "F")

# 6
# SAME
root_out3 = root_pcp_noncvx(X, lambda, mu, 5, verbose = TRUE)
#"Converged in 145 iterations."
norm(root_out3[[1]], "F")
norm(root_out3[[2]], "F")

# 5
# SAME
nan_nn_out = root_pcp_na_nonnegL(Xmissing, lambda, mu, verbose = TRUE)
norm(nan_nn_out[[1]], "F")
norm(nan_nn_out[[2]], "F")

# 4
# SAME
nan_nn_lod_out = root_pcp_na_nonnegL_lod(X, lambda, mu, matrix(0, nrow = nrow(X), ncol = ncol(X))
                                         , verbose = TRUE)
norm(nan_nn_lod_out[[1]], "F")
norm(nan_nn_lod_out[[2]], "F")

# 3
# SAME
n_out = root_pcp_nonnegL(X, lambda, mu, verbose = TRUE)
norm(n_out[[1]], "F")
norm(n_out[[2]], "F")

# 2
# SAME
nan_out = root_pcp_na(Xmissing, lambda, mu, verbose = TRUE)
norm(nan_out[[1]], "F")
norm(nan_out[[2]], "F")

# 1
# SAME
pcp_out = root_pcp(X, lambda, mu, verbose = TRUE)
norm(pcp_out[[1]], "F")
norm(pcp_out[[2]], "F")


