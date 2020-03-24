############################################################
# PCP w/ <LOD Penalty ######################################
# Original MATLAB by Jingkai ###############################
# 7/17/2019 ################################################
# 5 functions total ########################################
############################################################
# Updated: 2/19/2020, fixed matrix*matrix vs matrix*vector #
# Updated: 3/24/2020, removed negative pre-processing <LOD #
############################################################

library(matconv)
library(tidyverse)

############################################################
############################################################

# Prox L1 norm function, soft thresholding
# If Y < c (threshold), push to zero
prox_l1 <- function(Y, c) {
  
  myzero <- matrix(data = 0, ncol = ncol(Y), nrow = nrow(Y))
  X <- sign(Y) * pmax(abs(Y) - c, myzero, na.rm = TRUE)
  X
} 

############################################################
############################################################

# Prox nuclear norm function, L1 norm of the singular values
# This encourages matrix to be low rank by pushing SV to zero (sparse)
prox_nuclear <- function(Y,c) {
  
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
## this is kind of a MATLAB hack
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
loss_lod <- function(X, D, LOD) {
  # % D is the original data
  # % X is the new thing (L + S)
  # # LOD is the LOD
    X_lod <- (X - D)     * (D >= LOD) +
            ((X) - LOD)  * (D < LOD & ((X) > LOD)) +
              X          * (D < LOD & X < 0)
  
  l <- sum(X_lod^2) / 2
  # % L2 norm
  
  # % Any D_ij < 0 AND X_ij < LOD AND > 0 are treated as equal
  
  # % Minimize discrepancy for valid data
  # % Want to shrink negative things
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

pcp_lod <- function(D, lambda, mu, LOD) {
  
  if( any(D < matrix(0, nrow = nrow(D), ncol = ncol(D))) ) stop("Physical concentration measurements cannot be negative. GTFO.")
  
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
  MAX_ITER <- 100
  
  # Convergence Thresholds
  LOSS_THRESH <- 1e-5
  SAME_THRESH <- 1e-4
  
  if (is.vector(LOD)) {
    LOD = matrix(rep(LOD, times = nrow(D)), nrow = nrow(D), byrow = TRUE)
  }
  
  if( any(LOD < matrix(0, nrow = nrow(D), ncol = ncol(D))) ) stop("Physical Limits of Detection (LOD) cannot be negative. GTFO.")
  
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

      L2_new <- (L2_opt1 * (D >= LOD)) +
        (L2_opt2 * ((D < LOD) & ((L2 + S2) >= 0) & ((L2 + S2) <= LOD))) +
        (L2_opt3 * ((D < LOD) & ((L2 + S2) > LOD))) +
        (L2_opt4 * ((D < LOD) & ((L2 + S2) < 0)))

      S2_opt1 <- (mu*rho*D     + (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)
      S2_opt2 <- S1 + (Z3/rho)
      S2_opt3 <- (((mu*rho*LOD) + (((mu + rho)*Z3) - (mu*Z1) + ((mu + rho)*rho*S1) - (mu*rho*L1)))) / ((2*mu*rho) + (rho^2))
      S2_opt4 <- (               (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)

      S2 <- (S2_opt1 * (D >= LOD)) +
        (S2_opt2 * ((D < LOD) & ((L2 + S2) >= 0) & ((L2 + S2) <= LOD))) +
        (S2_opt3 * ((D < LOD) & ((L2 + S2) > LOD))) +
        (S2_opt4 * ((D < LOD) & ((L2 + S2) < 0)))
   
    L2 <- L2_new
    
    L3 <- pmax(L1 + Z2/rho, 0, na.rm = TRUE)
    # % Non-Negativity constraint!
    
    Z1 <- Z1 + rho*(L1 - L2)
    Z2 <- Z2 + rho*(L1 - L3)
    Z3 <- Z3 + rho*(S1 - S2)
    # % Z accumulate differnces between L and L and between S and S
    
    loss[i] <- nuclearL1 + 
      (lambda*sum(abs(S1))) +
      (mu*loss_lod((L2 + S2), D, LOD)) +
      sum(Z1*(L1 - L2)) +
      sum(Z2*(L1 - L3)) +
      sum(Z3*(S1 - S2)) +
      (rho/2 * (sum((L1-L2)^2) + sum((L1 - L3)^2) + sum((S1 - S2)^2)))
    # % The code block above takes LOD into account.
    
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

# x <- matrix(1:100, nrow = 20)
# pcp_lod(x, 1, 1, 0)
# 
# x[2,3] <- -0.0000001
# pcp_lod(x, 1, 1, 0)

# v <- 1:5
# if (is.vector(v)) {
#   LOD = matrix(rep(v, times = nrow(x)), nrow = nrow(x), byrow = TRUE)
# }

