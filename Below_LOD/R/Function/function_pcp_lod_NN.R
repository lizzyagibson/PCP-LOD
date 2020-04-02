############################################################
# PCP w/ <LOD Penalty ######################################
# Original MATLAB by Jingkai ###############################
# 7/17/2019 ################################################
# 5 functions total ########################################
############################################################
# Updated: 2/19/2020, fixed matrix*matrix vs matrix*vector #
# Updated: 3/34/2020, removed negative pre-processing <LOD #
############################################################

library(matconv)
library(tidyverse)

############################################################
############################################################

# Prox L1 norm function, soft thresholding
# if Y < c (threshold), push to zero
prox_l1 <- function(Y, c) {
  
  myzero <- matrix(data = 0, ncol = ncol(Y), nrow = nrow(Y))
  X <- sign(Y) * pmax(abs(Y) - c, myzero, na.rm = TRUE)
  X
} 

############################################################
############################################################

# Prox nuclear norm function, L1 norm of the singular values
# This encourages matrix to be low rank by pushing SV to zero (sparse)
prox_nuclear <- function(Y, c) {
  
  USV <- svd(Y)
  U <- USV$u
  S <- USV$d
  V <- USV$v
  
  myzero <- vector("numeric", length = length(S))
  S_new <- sign(S) * pmax(abs(S) - c, myzero, na.rm = TRUE)
  # Threshold the singular values, if SV < c, push it to zero
  
  X <- U %*% diag(S_new) %*% t(V)
  # % X is the truncation of the original
  # % Multiply the thresholded SVD components back together
  
  nuclearX  <- sum(abs(S_new))
  # This is the L1 norm of the truncated singular values
  # Goes into the loss function
  
  list(X = X, nuclearX = nuclearX)
}

############################################################
############################################################

# is same function for convergence criteria
# is the difference among matrices > noise threshold?
## if TRUE, keep iterating, if FALSE, end

# Compares L1, L2, L3 OR S1, S2
# Need ... for function to handle different number of inputs
# length(varargin) gives the number of function input arguments given in the call
# for L1, L2, L3, THREE comparisons, L1/L2, L1/L3, and L2/L3
# for S1, S2, one comparison
is_same <- function(SAME_THRESH, ...) {
  flag <- TRUE
  varargin <- list(...)
  if (length(varargin) == 2) {
    if (max(abs(varargin[[1]] - varargin[[2]])) > SAME_THRESH) {
      flag <- FALSE
    }
  }
  else if (length(varargin) == 3) {
    if ((max(abs(varargin[[1]] - varargin[[2]])) > SAME_THRESH) |
        (max(abs(varargin[[1]] - varargin[[3]])) > SAME_THRESH) |
        (max(abs(varargin[[2]] - varargin[[3]])) > SAME_THRESH)) {
      flag <- FALSE
    }
  }
  flag
}

############################################################
############################################################

# loss_lod function (only used in the loss function)
loss_lod_NN <- function(X, D, LOD) {
  # % D is the original data
  # % X is the new thing (L + S)
  # # LOD is the LOD
    X_lod <- (X - D)   * (D >= LOD) +
             (X - LOD) * ((D < LOD) & (X > LOD)) +
              X        * ((D < LOD) & (X < 0)) #+
             #(X - D)   * (D < LOD & (X > 0 && X <= LOD )) # or should it be zero
  
  l <- sum(X_lod^2) / 2
  # % L2 norm
  # % Any D_ij < 0 AND X_ij < LOD AND > 0 are treated as equal
  # % Minimize discrepancy for valid data
  # % Shrink negative things
  l
}

############################################################
############################################################

# % If the LOD threshold LOD = 0, solve the following ADMM splitting problem:
# %  min_{L1,L2,L3,S1,S2}
# %      ||L1||_* + lambda * ||S1||_1 + mu/2 * ||L2+S2-D||_F^2 + I_{L3>=0}
# %  s.t. L1 = L2
# %      L1 = L3
# %      S1 = S2.
# %
# % If LOD is not 0, replace ||L2+S2-D||_F^2 with LOD penalty.
# %
# % Below-LOD data input in D should be denoted as negative values, e.g. -1.

pcp_lod_NN <- function(D, lambda, mu, LOD) {
  
  if( any(D < matrix(0, nrow = nrow(D), ncol = ncol(D))) ) stop("Physical concentration measurements cannot be negative.")
  
  m <- nrow(D)
  n <- ncol(D)
  rho <- 1 # Augmented Lagrangian coefficient (rate)
  
  L1 <- matrix(0, m, n)
  L2 <- matrix(0, m, n)
  L3 <- matrix(0, m, n)
  
  S1 <- matrix(0, m, n)
  S2 <- matrix(0, m, n)
  
  Z1 <- matrix(0, m, n)
  Z2 <- matrix(0, m, n)
  Z3 <- matrix(0, m, n)
  
  # Max iteration
  MAX_ITER <- 1000
  
  # Convergence Thresholds
  LOSS_THRESH <- 1e-5
  SAME_THRESH <- 1e-4
  
  if (is.vector(LOD)) {
    t = matrix(TRUE, nrow = nrow(D), ncol = ncol(D))
    #LOD = t * LOD
    LOD = t(t(t) * delta10)
  }
  
  if( any(LOD < matrix(0, nrow = nrow(LOD), ncol = ncol(LOD))) ) stop("Physical Limits of Detection (LOD) cannot be negative.")
  
  loss <- vector("numeric", MAX_ITER)
  
  for (i in 1:MAX_ITER) {
    
    nuc <- prox_nuclear( ((L2 + L3 - (Z1 + Z2)/rho)/2), 1/2/rho)
    L1 <- nuc[[1]]
    nuclearL1 <- nuc[[2]] #nuclearX
    # % L, Z, S all start at zero, and change each iteration
    # % Prox_nuc is singular value thresholding
    # % L is low rank matrix
    
    S1 <- prox_l1((S2 - Z3/rho), lambda/rho)
    # % S is sparse matrix
    
      L2_opt1 <- (mu*rho*D     + (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)
      L2_opt2 <- L1 + Z1/rho
      L2_opt3 <- ((mu*rho*LOD + (((mu + rho)*Z1) - (mu*Z3) + ((mu + rho)*rho*L1) - (mu*rho*S1)))) / ((2*mu*rho) + (rho^2))
      L2_opt4 <- (               (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)

      L2_new <- (L2_opt1 * ( D >= LOD)) +
                (L2_opt2 * ((D < LOD) & (((L2 + S2) >= 0) & ((L2 + S2) <= LOD)))) +
                (L2_opt3 * ((D < LOD) & (((L2 + S2) > LOD)))) +
                (L2_opt4 * ((D < LOD) & (((L2 + S2) < 0))))

      S2_opt1 <- (mu*rho*D     + (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)
      S2_opt2 <- S1 + (Z3/rho)
      S2_opt3 <- (((mu*rho*LOD) + (((mu + rho)*Z3) - (mu*Z1) + ((mu + rho)*rho*S1) - (mu*rho*L1)))) / ((2*mu*rho) + (rho^2))
      S2_opt4 <- (               (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)

      S2 <- (S2_opt1 * ( D >= LOD)) +
            (S2_opt2 * ((D < LOD) & (((L2 + S2) >= 0) & ((L2 + S2) <= LOD)))) +
            (S2_opt3 * ((D < LOD) & (((L2 + S2) > LOD)))) +
            (S2_opt4 * ((D < LOD) & (((L2 + S2) < 0))))
   
    L2 <- L2_new
    
    L3 <- pmax(L1 + Z2/rho, 0, na.rm = TRUE)
    # % Non-Negativity constraint!
    
    Z1 <- Z1 + rho*(L1 - L2)
    Z2 <- Z2 + rho*(L1 - L3)
    Z3 <- Z3 + rho*(S1 - S2)
    # % Z accumulate differnces between L and L and between S and S
    
    loss[i] <- nuclearL1 + 
      (lambda*sum(abs(S1))) +
      (mu*loss_lod_NN((L2 + S2), D, LOD)) +
      sum(Z1*(L1 - L2)) +
      sum(Z2*(L1 - L3)) +
      sum(Z3*(S1 - S2)) +
      (rho/2 * (sum((L1-L2)^2) + sum((L1 - L3)^2) + sum((S1 - S2)^2)))
    # % The code block above takes LOD into account.
    
    print(str_c(i, " Obj: ", loss[i]))
    
    if ((i != 1) && 
        (abs(loss[i-1] - loss[i]) < LOSS_THRESH) && 
        is_same(SAME_THRESH, L1, L2, L3) &&
        is_same(SAME_THRESH, S1, S2)) {
      break} # % Convergence criteria!
  }
  
  L <- L3 # (L1 + L2 + L3) / 3
  S <- S1 #(S1 + S2) / 2
  list(L = L, S = S, loss = loss)
}

# Compare with original <LOD version
# x <- matrix(1:100, nrow = 20)
# nn_out <- pcp_lod_NN(x, 1, 1, 0)
# reg_out <- pcp_lod(x, 1, 1, 0)
# 
# colSums(nn_out$L)
# colSums(reg_out$L)

# # Read air pollution data
# library(R.matlab)
# # mixture <- readMat(here::here("Data/mixtures_data.mat"))
#  
# mix <- as.data.frame(mixture) %>% as_tibble() %>%
#   select(Al, As, Ba, bc, Br, Ca, Cl,
#          Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
#          Ti,  V, Zn) %>%
#   drop_na(.)
# 
# m <- nrow(mix)
# n <- ncol(mix)
# lambda_mix = 1/sqrt(m)
#  
# # m_nn_out <- pcp_lod_NN(mix, lambda_mix, 10, 0)
# # m_reg_out <- pcp_lod(mix, lambda_mix, 10, 0)
# 
# # colSums(m_nn_out$L)
# # colSums(m_reg_out$L)
# 
# # Increasing
# mix <- as_tibble(mix)
# 
# mix_data_lod_10 <- mix %>% 
#   mutate_all(~ifelse(. < quantile(., probs = .10), 0, .)) %>% as.matrix()
# 
# mix_data_lod_20 <- mix %>% 
#   mutate_all(~ifelse(. < quantile(., probs = .20), 0, .)) %>% as.matrix()
# 
# mix_data_lod_30 <- mix %>% 
#   mutate_all(~ifelse(. < quantile(., probs = .30), 0, .)) %>% as.matrix()
# 
# mix_data_lod_40 <- mix %>% 
#   mutate_all(~ifelse(. < quantile(., probs = .40), 0, .)) %>% as.matrix()
# 
# mix_data_lod_50 <- mix %>% 
#   mutate_all(~ifelse(. < quantile(., probs = .50), 0, .)) %>% as.matrix()
# 
# delta10 <- mix %>% 
#   summarise_all(quantile, probs = .10) %>% as_vector()
# 
# delta20 <- mix %>% 
#   summarise_all(quantile, probs = .20) %>% as_vector()
# 
# delta30 <- mix %>% 
#   summarise_all(quantile, probs = .30) %>% as_vector()
# 
# delta40 <- mix %>% 
#   summarise_all(quantile, probs = .40) %>% as_vector()
# 
# delta50 <- mix %>% 
#   summarise_all(quantile, probs = .50) %>% as_vector()
# 
# mix <- as.matrix(mix)
# results_0  <- pcp_lod(mix,             lambda_mix, 10, 0)
# results_10 <- pcp_lod(mix_data_lod_10, lambda_mix, 10, delta10)
# results_20 <- pcp_lod(mix_data_lod_20, lambda_mix, 10, delta20)
# results_30 <- pcp_lod(mix_data_lod_30, lambda_mix, 10, delta30)
# results_40 <- pcp_lod(mix_data_lod_40, lambda_mix, 10, delta40)
# results_50 <- pcp_lod(mix_data_lod_50, lambda_mix, 10, delta50)
# 
# nn_results_0  <- pcp_lod_NN(mix, lambda_mix, 10, 0)
# nn_results_10 <- pcp_lod_NN(mix_data_lod_10, lambda_mix, 10, delta10)
# nn_results_20 <- pcp_lod_NN(mix_data_lod_20, lambda_mix, 10, delta20)
# nn_results_30 <- pcp_lod_NN(mix_data_lod_30, lambda_mix, 10, delta30)
# nn_results_40 <- pcp_lod_NN(mix_data_lod_40, lambda_mix, 10, delta40)
# nn_results_50 <- pcp_lod_NN(mix_data_lod_50, lambda_mix, 10, delta50)
# 
# all.equal(colSums(results_0$L),  colSums(nn_results_0$L))
# all.equal(colSums(results_10$L), colSums(nn_results_10$L))
# all.equal(colSums(results_20$L), colSums(nn_results_20$L))
# all.equal(colSums(results_30$L), colSums(nn_results_30$L))
# all.equal(colSums(results_40$L), colSums(nn_results_40$L))
# all.equal(colSums(results_50$L), colSums(nn_results_50$L))
# 
# # Diff negative pre-processing gives same results
# # mix_data_lod_50_100 <- as_tibble(mix) %>%
# #   mutate_all(~ifelse(. <= quantile(., probs = .50), -100, .)) %>% as.matrix()
# # results_50_1 <-   pcp_lod(mix_data_lod_50, lambda_mix, 10, delta50)
# # results_50_100 <- pcp_lod(mix_data_lod_50_100, lambda_mix, 10, delta50)
# # all.equal(colSums(results_50_1$L),  colSums(results_50_100$L))


