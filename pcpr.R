
# PCP-LOD function and subfunctions ####

## helper functions ####

# Prox nuclear
#
# \code{prox_nuclear} implements the proximal gradient method for the nuclear norm.
# The nuclear norm is equivalent to the L2 norm of the singular values of the matrix.
# This singular value thresholding encourages the \code{L} matrix to be low rank.
# The proximal gradient method is used to solve non-differentiable convex optimization problems.
# This is used in \code{stablePCP} and \code{rootPCP}.
#
# @param Y The \code{L} matrix.
# @param c The amount by which the prox nuclear method penalizes \code{Y}.
#
# @return \code{prox_nuclear} returns two objects.
# \code{X} is the thresholded \code{L} matrix.
# \code{nuclearX} is the sum of the absolute values of the thresholded singular values.
# This goes into the objective function.
prox_nuclear <- function(Y,c) {
  
  USV <- svd(Y)
  U <- USV$u
  S <- USV$d
  V <- USV$v
  
  S_new <- sign(S) * pmax(abs(S) - c, 0)
  X <- U %*% diag(S_new) %*% t(V)
  nuclearX  <- sum(abs(S_new))
  
  return(list(X = X, nuclearX = nuclearX))
}

# Prox L1
#
# \code{prox_l1} implements the proximal gradient method for the L1 norm.
# This soft thresholding encourages the \code{S} matrix to be sparse.
# The proximal gradient method is used to solve non-differentiable convex optimization problems.
# This is used in \code{stablePCP} and \code{rootPCP}.
#
# @param Y The \code{S} matrix.
# @param c The amount by which the prox L1 method penalizes \code{Y}.
#
# @return The thresholded \code{S} matrix.
#
prox_l1 <- function(Y, c) {
  X <- sign(Y) * pmax(abs(Y) - c, 0)
  return(X)
}

# Prox Frobenius
#
# \code{prox_fro} implements the proximal gradient method for the Frobenius norm.
# This thresholding minimizes the square root of the sum of squared error.
# The proximal gradient method is used to solve non-differentiable convex optimization problems.
# This is only used in \code{rootPCP} because the squared Frobenius error in \code{stablePCP} is differentiable.
#
# @param X The error matrix, \code{D-L-S}
# @param c The amount by which the prox Frobenius method penalizes \code{X}.
#
# @return The thresholded error matrix.
#
prox_fro <- function(X,c) {
  
  n = norm(X,'F')
  
  if (n <= c) {Y = matrix(0, nrow = nrow(X), ncol = ncol(X))}
  else {Y = (1 - c/n) * X}
  
  return(Y)
}

# Loss for stablePCP-LOD
#
# \code{loss_lod} includes the LOD-specific penalty terms to compute the loss on the squared error term \code{L+S-D}.
#
# @param X The predicted value of \code{L + S} at the current iteration.
# @param D The original dataset.
# @param LOD The LOD. May be a scalar, vector (\code{length(LOD) = ncol(D)}), or matrix (\code{dim(LOD) == dim(D)}).
#
# @return Scalar value used to calculate loss in objective function.
#
loss_lod <- function(X, D, LOD) {
  X_lod <- (X - D)   * (D >= 0) +
    (X - LOD)  * ((D < 0) & (X > LOD)) +
    X        * ((D < 0) & (X < 0))
  
  sum(X_lod^2) / 2
}

# Project rank
# Non-convex replacement for nuclear norm
#
proj_rank_r = function(Y, r) {
  
  if (ncol(Y) < r) stop("r > matrix rank")
  if (ncol(Y) == r) {return(Y)}
  
  USV <- svd(Y)
  U <- USV$u
  S <- USV$d
  V <- USV$v
  
  s = S
  
  # used to read:
  #if (length(S)-1 == r) {
  #  s[r]  = 0
  #} else {
  #  s[(r+1):length(s)] = 0
  #}
  # now reads:
  s[(r+1):length(s)] = 0
  
  S_new  = diag(s)
  X = U %*% S_new %*% t(V)
  return(X)
}

# HT
# Hard-thresholding for the sparse matrix
#
HT = function(Y, c) {
  X = Y
  X[abs(X) < c] = 0
  return(X)
}

# PCP-LOD ####

# Nonconvex nonnegative squareroot PCP function with missing values (NA) & LOD penalty
#
# \code{root_pcp_noncvx_nonnegL_na_lod} implements nonconvex \code{rootPCP} with a non-negativity constraint on
# the \code{L} solution matrix and LOD-specific penalties. \cr \cr
# It solved the following ADMM splitting problem: \cr \cr
# min(L,S) \cr
# 1_{rank(L) <= r} + lambda * ||S||_1 + mu * ( sum_(ij observed) f((L+S)_ij, D_ij) )^0.5 + I[L>=0] \cr \cr
# This is first transformed to the problem: \cr \cr
# min(L1,L2,S,Z) \cr
# 1_{rank(L1) <= r} + lambda * ||S1||_1 + mu * ( sum_(ij observed) f(Z_ij,D_ij) )^0.5 + I[L3>=0] \cr \cr
# s.t. L1 = L2; Z = P_obs[ L2 + S2]; L1 = L3 \cr \cr
# The algorithm conducts ADMM splitting as (L1, S1, Z), (L2, S2, L3). \cr \cr
# This version allows for missing values. \cr \cr
# Use NA for missing entries in D. \cr \cr
# Assume that the true L >= 0 & observations in D >= 0 \cr \cr
# Use -1 for <LOD \cr \cr
# LOD penalty: \cr \cr
# f(x,y) = (x-y)^2 if y>0 \cr
#        = (x-LOD)^2 if y=-1, x>LOD \cr
#        = x^2 if y=-1, x<0 \cr
#        = 0 otherwise
#
# @param D The original dataset.
# @param lambda The \code{lambda} parameter penalizes the proximal L1 gradient on the \code{S} matrix.
# @param mu The \code{mu} parameter penalizes the error term.
# @param r The \code{r} parameter specifies the desired rank.
# @param LOD The LOD (limit of detection) may be a scalar, vector (\code{length(LOD) = ncol(D)}), or matrix (\code{dim(LOD) == dim(D)}).
# @param verbose A logical indicating if you would like information on the number of iterations required to reach convergence printed. Optional, and by default \code{verbose = FALSE}.
# @param MAX_ITER optional parameter to set the maximum number of iterations, default = 20,000
#
# @return Returns two solution matrices, the low rank \code{L} matrix and the sparse \code{S} matrix.
#
# @export
root_pcp_noncvx_nonnegL_na_lod <- function(D, lambda, mu, r, LOD, verbose = FALSE, MAX_ITER = 10000) {
  
  if (any(class(LOD) == "list")) {
    LOD <- unlist(LOD)
  }
  
  n = nrow(D)
  p = ncol(D)
  rho = 0.1; # Augmented Lagrangian parameter
  
  L1 <- matrix(0, n, p)
  L2 <- matrix(0, n, p)
  L3 <- matrix(0, n, p)
  
  S1 <- matrix(0, n, p)
  S2 <- matrix(0, n, p)
  
  Z  <- matrix(0, n, p)
  Y1 <- matrix(0, n, p)
  Y2 <- matrix(0, n, p)
  Y3 <- matrix(0, n, p)
  Y4 <- matrix(0, n, p)
  
  # mask: support of observation of D
  
  mask_above_lod = D >= 0
  mask_above_lod[is.na(mask_above_lod)] = 0
  
  mask_below_lod = D < 0
  mask_below_lod[is.na(mask_below_lod)] = 0
  
  mask_obs = !is.na(D)
  D[!mask_obs] = -2
  
  EPS_ABS = 1e-6
  EPS_REL = 1e-6
  
  flag_converge = 0
  
  if (is.vector(LOD) & length(LOD) != 1) {
    LOD = kronecker(matrix(1,n),t(LOD))
  } # This converts a vector LOD to a matrix, so that it multiplies correctly
  
  #% ADMM-splitting iterations
  for (i in 1:MAX_ITER) {
    
    #% Store previous values of L2,S<Z
    L2_old = L2
    S2_old = S2
    L3_old = L3
    
    #% Update 1st primal variable (L1,S1,Z)
    L1 = proj_rank_r( (L2+L3-Y1/rho-Y4/rho)/2, r)
    
    S1 = prox_l1( S2-Y2/rho, lambda/rho )
    
    temp = L2+S2-Y3/rho
    # Z_unobs = (1-mask_obs) * temp
    # Z_obs_below_LOD1 = (mask_below_lod & (temp>=0) & (temp<=LOD)) * temp
    # temp2 = mask_obs * (1-(mask_below_lod & (temp>=0) & (temp<=LOD))) * temp -
    #         (mask_above_lod * D) - (LOD * mask_below_lod * (temp>=LOD))
    # Z = prox_fro( temp2, mu/rho ) + (mask_above_lod * D) +
    #    (LOD * mask_below_lod * (temp>=LOD)) +
    #    Z_unobs + Z_obs_below_LOD1
    
    temp_D = D*mask_above_lod +
      temp*(mask_below_lod & (temp>=0) & (temp<=LOD)) +
      LOD*(mask_below_lod & (temp>LOD))
    Z = prox_fro( temp - temp_D, mu/rho ) + temp_D
    
    L3 = pmax(L1+Y4/rho,0)
    
    #% Update 2nd primal variable (L2,S2)
    L2_obs = mask_obs * (1/3*( 2*L1-S1+Z + (2*Y1-Y2+Y3)/rho ))
    L2_unobs = (1-mask_obs)*(L1+Y1/rho)
    L2 = L2_obs+L2_unobs
    
    S2_obs = mask_obs * (1/3*( 2*S1-L1+Z + (2*Y2-Y1+Y3)/rho ))
    S2_unobs = (1-mask_obs) * (S1+Y2/rho)
    S2 = S2_obs+S2_unobs
    
    #% Update dual variable (Y1,Y2)
    Y1 = Y1 + rho*(L1-L2)
    Y2 = Y2 + rho*(S1-S2)
    Y3 = Y3 + rho*mask_obs*(Z-(L2 + S2))
    Y4 = Y4 + rho*(L1-L3)
    
    #%  Calculate primal & dual residuals; Update rho
    res_primal = sqrt(norm(L1-L2,'F')^2 + norm(S1-S2,'F')^2 + norm(mask_obs*(Z-L2-S2),'F')^2 + norm(L1-L3,'F')^2)
    res_dual = rho * sqrt( norm(L2+L3-L2_old-L3_old,'F')^2 + norm(S2-S2_old,'F')^2 +
                             norm(mask_obs*(L2-L2_old+S2-S2_old),'F')^2 )
    
    if (res_primal > 10 * res_dual) {
      rho = rho * 2
    } else if (res_dual > 10 * res_primal) {
      rho = rho / 2}
    
    #% Check stopping criteria
    thresh_primal = EPS_ABS * sqrt(4*n*p) + EPS_REL *
      max(sqrt( 2*norm(L1,'F')^2 + norm(S1,'F')^2 + norm(Z,'F')^2 ),
          sqrt( norm(L2,'F')^2 + norm(S2,'F')^2 + norm(mask_obs*(L2+S2),'F')^2 + norm(L3,'F')^2))
    thresh_dual = EPS_ABS * sqrt(3*n*p) + EPS_REL *
      sqrt( norm(Y1+Y4,'F')^2 + norm(Y2,'F')^2 + norm(Y3,'F')^2  )
    
    final_iter = i
    if (res_primal < thresh_primal && res_dual < thresh_dual) {
      flag_converge = 1
      if (verbose) print(paste0('Converged in ', i,' iterations.'))
      break}
    
  }
  
  L = (L1+L2+L3)/3
  S = (S1+S2)/2
  
  if (flag_converge == 0 & verbose) print('Did not converge.')
  L[L < 0] <- 0
  return(list(L = L, S = S, final_iter = final_iter))
}

