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
