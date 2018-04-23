#MATLAB CODE

#function w = soft_thresholding( v, lambda )
#w = sign(v) .* max( abs(v) - lambda, 0 ); 

#% .* is element wise multiplication
#% sign of v, preserve the positive sign of the singular values
#% shrink the value in v toward zero by lambda
#% take the higher of zero of value minus lambda

#R code

soft_thresholding <- function(v, lambda) {
  myzero <- matrix(data = 0, ncol = ncol(v), nrow = nrow(v))
  w <- sign(v) * pmax(abs(v) - lambda, myzero)
  w
} 

soft_thresholding_diag <- function(v, lambda) {
  myzero <- vector("numeric", length = length(v))
  w <- sign(v) * pmax(abs(v) - lambda, myzero)
  w
} 
