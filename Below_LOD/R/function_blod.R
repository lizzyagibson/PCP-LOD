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
  
    myzero <- matrix(data = 0, ncol = ncol(v), nrow = nrow(v))
    X <- sign(v) * pmax(abs(v) - lambda, myzero)
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
  S_new <- sign(S) * pmax(abs(S) - c, myzero)
  # Threshold the singular values, if SV < c, push it to zero
  
  X <- U %*% S_new %*% t(V)
  
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
  
  X_lod <- (X - D) * (D >= 0) +
    (X - Delta) * (D < 0 & X > Delta) +
    X * (D < 0 & X < 0)
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





