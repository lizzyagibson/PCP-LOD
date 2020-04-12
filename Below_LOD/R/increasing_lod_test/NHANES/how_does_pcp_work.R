mix_50 <- mix_data_lod_50
mix_50[,1] <- -1

check_pcp <- pcp_lod_j(mix_50, lambda = lambda_mix, mu = mu_mix, LOD = delta50)

head(L_lod50)[,1:10]
head(L_lod50_sqrt2)[,1:10]
head(mix_50)[,1:10]
head(check_pcp$change[,,100])

check_pcp$change[3,10,]

pcp_lod_j <- function(D, lambda, mu, LOD) {
  
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
    empty = matrix(1, nrow = nrow(D), ncol = ncol(D))
    LOD = t(t(empty) * LOD)
  }
  
  loss <- vector("numeric", MAX_ITER)
  change <- array(NA, c(m, n, MAX_ITER))
  
  for (i in 1:MAX_ITER) {
    
    nuc <- prox_nuclear( ((L2 + L3 - (Z1 + Z2)/rho)/2), 1/2/rho)
    L1 <- nuc[[1]]
    nuclearL1 <- nuc[[2]] #nuclearX
    
    S1 <- prox_l1((S2 - Z3/rho), lambda/rho)
    
    L2_opt1 <- (mu*rho*D     + (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)
    L2_opt2 <- L1 + Z1/rho
    L2_opt3 <- ((mu*rho*LOD + (((mu + rho)*Z1) - (mu*Z3) + ((mu + rho)*rho*L1) - (mu*rho*S1)))) / ((2*mu*rho) + (rho^2))
    L2_opt4 <- (               (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)
    
    L2_new <- (L2_opt1 * (D >= 0)) +
      (L2_opt2 * ((D < 0) & (((L2 + S2) >= 0) & ((L2 + S2) <= LOD)))) +
      (L2_opt3 * ((D < 0) & (((L2 + S2) > LOD)))) +
      (L2_opt4 * ((D < 0) & (((L2 + S2) < 0))))
    
    S2_opt1 <- (mu*rho*D     + (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)
    S2_opt2 <- S1 + (Z3/rho)
    S2_opt3 <- (((mu*rho*LOD) + (((mu + rho)*Z3) - (mu*Z1) + ((mu + rho)*rho*S1) - (mu*rho*L1)))) / ((2*mu*rho) + (rho^2))
    S2_opt4 <- (               (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)
    
    S2 <- (S2_opt1 * (D >= 0)) +
      (S2_opt2 * (((D < 0) & ((L2 + S2) >= 0) & ((L2 + S2) <= LOD)))) +
      (S2_opt3 * (((D < 0) & ((L2 + S2) > LOD)))) +
      (S2_opt4 * (((D < 0) & ((L2 + S2) < 0))))
    
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
    
    change[,,i] <- L3
    print(str_c(i, " Obj: ", loss[i]))
    
    if ((i != 1) && 
        (abs(loss[i-1] - loss[i]) < LOSS_THRESH) && 
        is_same(SAME_THRESH, L1, L2, L3) &&
        is_same(SAME_THRESH, S1, S2)) {
      break} # % Convergence criteria!
  }
  
  L <- L3 # (L1 + L2 + L3) / 3
  S <- S1 #(S1 + S2) / 2
  list(L = L, S = S, loss = loss, change = change)
}
