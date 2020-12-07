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
# X[1:5, 1:5]

m <- nrow(X)
n <- ncol(X)

## Run 3 models
#tic()
lod_out = pcp_lod(X, 1/sqrt(m), 10, 0)
#"Converged in 69 iterations."
#toc()
norm(lod_out[[1]], "F")
norm(lod_out[[2]], "F")

tic()
root_out = root_pcp(X, 1/sqrt(m), 10)
#"Converged in 837 iterations."
toc()
# 14.049 sec elapsed
norm(root_out[[1]], "F")
norm(root_out[[2]], "F")

#tic()
root_nn = root_pcp_nonnegL(X, 1/sqrt(m), 10)
#"Converged in 646 iterations."
#toc()
norm(root_nn[[1]], "F")
norm(root_nn[[2]], "F")

tic()
nan_out = root_pcp_with_nan(X, 1/sqrt(m), 10)
#"Converged in 837 iterations."
toc()
# 16.649 sec elapsed
norm(nan_out[[1]], "F")
norm(nan_out[[2]], "F")

## With missing values!
Xmissing = X
Xmissing[1:1000,1:5] = NA

tic()
nan_out = root_pcp_with_na(Xmissing, 1/sqrt(m), 10)
#"Converged in 915 iterations."
toc()
# 18.347 sec elapsed
norm(nan_out[[1]], "F")
norm(nan_out[[2]], "F")
