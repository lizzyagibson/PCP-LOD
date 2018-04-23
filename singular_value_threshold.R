#MATLAB code

#function [N,v,sv] = singular_value_threshold( M, lambda, B, sv )
#% B specified as [] in example
#% singular value thresholding
#% solves two problems

#[U,S,V] = svd(M);
#% SVD function in R
#% D in R is S in MATLAB
#% U is U and V is V

#N = U * soft_thresholding(S,lambda) * V'; 
#%prime is transpose
#%S is singular values

#end
#v = sum(diag(soft_thresholding(S,lambda))); 
#% diag converts diagonal matrix into straight vector
#% takes singular value diagonal matrix and makes vector of singular values

#R code

source("./PCP/soft_thresholding.R")

singular_value_threshold <- function(M, lambda) {
  
  USV <- svd(M)
  U <- USV$u
  S <- USV$d
  V <- USV$v
    
  N <- U %*% diag(soft_thresholding_diag(S, lambda)) %*% t(V)
  
  v  <- sum(soft_thresholding_diag(S, lambda))
  svt <- list(N = N, v = v) #is this N the same as L??
  svt
}
