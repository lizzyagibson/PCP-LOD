library(R.matlab)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(pcpr)

#####
##### Data
#####

mixture <- readMat(here::here("Data/mixtures_data.mat"))

X <- as.data.frame(mixture) %>% as_tibble() %>% 
  select(Al, As, Ba, bc, Br, Ca, Cl,
         Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
         Ti,  V, Zn) %>% 
  drop_na() %>% as.matrix()

X_std <- apply(X, 2, function(x) x/sd(x))

m = nrow(X)
n = ncol(X)

lambda = 1/sqrt(m)
mu = sqrt(n/(2*log(m*n)))

system.time(svd(X))
system.time(propack.svd(X))
system.time(fast.svd(X))

#####
##### Rho = 1
#####
# Run model reg
system.time(pcp_lod(X, lambda, mu, 0))
# % Elapsed time is 42.797 seconds.
# % Converges in 747
# PROPACK SVD converges in 28.444 seconds.
# Fast SVD in 26.960...

# % % Run model standardized
system.time(pcp_lod(X_std, lambda, mu, 0))
# % Elapsed time is 88.398 seconds.
# % Converges in 1905
# PROPACK SVD converges in 79.193
# FAST SVD converges in 69.711

#####
##### Rho = 0.025
#####

# % % Run model reg
system.time(pcp_lod(X, lambda, mu, 0))
# % Elapsed time is 15.941 seconds.
# % Converges in 383

# % % Run model standardized
system.time(pcp_lod(X_std, lambda, mu, 0))
# % Elapsed time is 482.900 seconds.
# % DOES NOT CONVERGE in 10,000

#####
##### Rho changes each iter
#####

# % % Run model reg
system.time(pcp_lod(X, lambda, mu, 0))
# % Elapsed time is 2.582 seconds.
# % Converges in 66

# % % Run model standardized
system.time(pcp_lod(X_std, lambda, mu, 0))
# % Elapsed time is 6.522 seconds.
# % Converges in 185

#####
##### Rho changes each iter + slow SVD
#####

# % % Run model reg
system.time(pcp_lod(X, lambda, mu, 0))
# % Elapsed time is 2.578 seconds.
# % Converges in 66

# % % Run model standardized
system.time(pcp_lod(X_std, lambda, mu, 0))
# % Elapsed time is 6.789 seconds.
# % Converges in 185
