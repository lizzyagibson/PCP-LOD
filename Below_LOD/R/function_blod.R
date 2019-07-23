############################################################
# PCP w/ <LOD Penalty ######################################
# Original MATLAB by Jingkai ###############################
# 7/17/2019 ################################################
# 5 functions total ########################################
############################################################

library(matconv)
library(tidyverse)

############################################################
############################################################

# Prox L1 norm function, soft thresholding
# if Y < c, push to zero
prox_l1 <- function(Y, c) {
  
    myzero <- matrix(data = 0, ncol = ncol(Y), nrow = nrow(Y))
    X <- sign(Y) * pmax(abs(Y) - c, myzero, na.rm = TRUE)
    X
    } 

############################################################
############################################################

# Prox nuclear norm function, L1 norm of the singular values
# This encourages matrix to be low rank by pushing SV to zero
prox_nuclear <- function(Y,c) {
  
  USV <- svd(Y)
  U <- USV$u
  S <- USV$d
  V <- USV$v

  myzero <- vector("numeric", length = length(S))
  S_new <- sign(S) * pmax(abs(S) - c, myzero, na.rm = TRUE)
  # Threshold the singular values, if SV < c, push it to zero
  
  X <- U %*% diag(S_new) %*% t(V)
  
  nuclearX  <- sum(abs(S_new))
  # This is the L1 of the truncated singular values
 
  list(X = X, nuclearX = nuclearX)
  }

############################################################
############################################################

# is same function for convergence criteria
# is the difference among matrices > noise threshold?
## if TRUE, keep iterating, if FALSE, end
is_same <- function(SAME_THRESH, ...) {
  flag <- TRUE
  varargin <- list(...)

  for (i in 1:(length(as.list(match.call())) - 3)) {
    if (max(abs(varargin[[i]] - varargin[[i + 1]])) > SAME_THRESH) {
        flag <- FALSE
      }
    }
  flag
}

############################################################
############################################################

# loss_lod function

loss_lod <- function(X, D, Delta) {
  
  X_lod <- (X - D)     * (D >= 0) +
           (X - Delta) * (D < 0 & X > Delta) +
            X          * (D < 0 & X < 0)
  # Element-wise boolean operator
  # If D_ij >= 0, then X_lod = X_ij - D_ij
  # If D_ij < 0 (This means < LOD) AND X_ij > Delta, then X_lod = X_ij - Delta
  # If D_ij < 0 AND X_ij < 0, then X_lod = X_ij
  l <- sum(X_lod^2) / 2
  # L2 norm
  # Want to minimize this loss (e.g., discrepancy from valid data), shrink negative things
  l
  }

############################################################
############################################################

# % If the LOD threshold Delta = 0, solve the following ADMM splitting problem:
#   % min_{L1,L2,L3,S1,S2}
# %      ||L1||_* + lambda * ||S1||_1 + mu/2 * ||L2+S2-D||_F^2 + I_{L3>=0}
# % s.t. L1 = L2
# %      L1 = L3
# %      S1 = S2.
# %
# % If Delta is not 0, replace ||L2+S2-D||_F^2 with LOD penalty.
# %
# % Below-LOD data input in D should be denoted as negative values, e.g. -1.

pcp_lod <- function(D, lambda, mu, Delta) {
  
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
  
  MAX_ITER <- 100
  LOSS_THRESH <- 1e-5
  SAME_THRESH <- 1e-4
  
  loss <- vector("numeric", MAX_ITER)
  
  for (i in 1:MAX_ITER) {
    
    nuc <- prox_nuclear( ((L2 + L3 - (Z1 + Z2)/rho)/2), 1/2/rho)
    L1 <- nuc[[1]]
    nuclearL1 <- nuc[[2]]
    
    S1 <- prox_l1( (S2 - Z3/rho), lambda/rho)
    
    L2_opt1 <- (mu*rho*D     + (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)
    L2_opt2 <- L1 + Z1/rho
    L2_opt3 <- (mu*rho*Delta + (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)
    L2_opt4 <- (               (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)
    
    L2_new <- (L2_opt1 * (D >= 0)) +
              (L2_opt2 * ((D < 0) & ((L2 + S2) >= 0) & ((L2 + S2) <= Delta))) +
              (L2_opt3 * ((D < 0) & ((L2 + S2) > Delta))) +
              (L2_opt4 * ((D < 0) & ((L2 + S2) < 0)))
    
    S2_opt1 <- (mu*rho*D     + (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)
    S2_opt2 <- S1 + (Z3/rho)
    S2_opt3 <- (mu*rho*Delta + (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)
    S2_opt4 <- (               (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)
    
    S2 <- (S2_opt1 * (D >= 0)) +
          (S2_opt2 * ((D < 0) & ((L2 + S2) >= 0) & ((L2 + S2) <= Delta))) +
          (S2_opt3 * ((D < 0) & ((L2 + S2) > Delta))) +
          (S2_opt4 * ((D < 0) & ((L2 + S2) < 0)))
  
    L2 <- L2_new
    
    L3 <- pmax(L1 + Z2/rho, 0, na.rm = TRUE)
    
    Z1 <- Z1 + rho*(L1 - L2)
    Z2 <- Z2 + rho*(L1 - L3)
    Z3 <- Z3 + rho*(S1 - S2)
    
    loss[i] <- nuclearL1 + 
      (lambda*sum(abs(S1))) +
      (mu*loss_lod((L2 + S2), D, Delta)) +
      sum(Z1*(L1 - L2)) +
      sum(Z2*(L1 - L3)) +
      sum(Z3*(S1 - S2)) +
      (rho/2 * (sum((L1-L2)^2) + sum((L1 - L3)^2) + sum((S1 - S2)^2)))
    
    if ((i != 1) && 
        (abs(loss[i-1] - loss[i]) < LOSS_THRESH) && 
        is_same(SAME_THRESH, L1, L2, L3) &&
        is_same(SAME_THRESH, S1, S2)) {
      break}
  }
  
  L <- (L1 + L2 + L3) / 3
  S <- (S1 + S2) / 2
  list(L = L, S = S, loss = loss)
}







