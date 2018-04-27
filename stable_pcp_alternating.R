#R code

library(stringr)

source("./PCP/singular_value_threshold.R")

stable_pcp_alternating <- function(D, lambda, mu) {
  
  m <- nrow(D)
  n <- ncol(D)
  
  S <- matrix(0, nrow = m, ncol = n)
  L <- matrix(0, nrow = m, ncol = n)
  
  iter <- 0
  MAX_ITER <- 20
  done <- FALSE
  
  while (!done) {
    
    iter <- iter + 1
    
    svt <- singular_value_threshold((D - S), 1/mu)
    L <- svt[[1]] #svt$N #Are N and L equivalent?
    v <- svt[[2]]
    
    S <- soft_thresholding((D - L), lambda/mu)
    
    obj <- v + lambda * sum(abs(S)) + (mu/2) * norm((D - L - S), type = "F")^2
    
    print(str_c(iter, " Obj: ", obj))
    
    if (iter >= MAX_ITER) {done <- TRUE}
    
  }
  list(L = L, S = S, Lambda = lambda, Mu = mu)
}

#More iterations, diff values for lambda will get lower rank L matrix