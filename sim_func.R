sim_data <- function(sim_seed, nrow, ncol, rank, sigma, add_sparse=FALSE) {
  
  # gaussian noise = Z:
  set.seed(994 + sim_seed)
  Z <- rnorm(nrow, ncol) * sigma
  
  # low rank matrix = L:
  set.seed(1996 + sim_seed)
  U <- rand(nrow, rank)
  set.seed(2998 + sim_seed)
  V <- rand(rank, ncol)
  L <- U%*%V # ground truth low rank matrix.
  
  # intermediate step to help define sparse matrix S below:
  D <- L + Z
  
  # sparse matrix = S:
  if (add_sparse) {
    set.seed(sim_seed)
    S <- -D * ((D < 0) * 1) + (rand(nrow,ncol)<0.03) * rand(nrow,ncol)*1 # Jingkai's method of simulating sparse noise (1st term ensures nonnegativity, second term adds some sparse noise to random entries)
    #S <- (-D + rand(nrow, ncol)) * ((D < 0) * 1) # Lawrence's method of simulating sparse noise (use sparse to ensure matrix is non-neg + add extra noise on top of those same entries)
    #S <- rsparsematrix(nrow = nrow, ncol = ncol, density = .07, rand.x = rand) # for when there is no gaussian (Z) noise added. buggy at the moment throwing fatal error.
  } else {
    S <- zeros(nrow, ncol)
  }
  
  # final simulated data matrix = M:
  M <- S + D
  M[M < 0] <- 0 # non-negative
  ret <- list(M = M, L = L, S = S, Z = Z)
  ret
}