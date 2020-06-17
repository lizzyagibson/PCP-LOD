############################################################
# PCP w/ <LOD Penalty ######################################
# Troubleshooting ##########################################
# 5 functions total ########################################
############################################################

library(matconv)
library(tidyverse)
library(R.matlab)

# Read air pollution data
mixture <- readMat("./Data/mixtures_data.mat")

mixture_data <- as.data.frame(mixture) %>% as_tibble() %>% 
  select(Al, As, Ba, bc, Br, Ca, Cl,
         Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
         Ti,  V, Zn) %>% 
  drop_na()

x <- as.matrix(mixture_data)

# Creat 10 x 10 matrix
dat <- matrix(1:100, nrow = 10, byrow = TRUE)

# Create 11 x 10 matrix
dat_11 <- rbind(dat, 101:110)

# 50% <LOD
mix_data_lod_50 <- mixture_data %>% 
  mutate_all(~ifelse(. <= quantile(., probs = .50), -1, .)) %>% as.matrix()

delta50 <- mixture_data %>% 
  summarise_all(quantile, probs = .50) %>% as_vector()

############################################################
############################################################

# Prox L1 norm function, soft thresholding
# if Y < c (threshold), push to zero
prox_l1 <- function(Y, c) {
  
    myzero <- matrix(data = 0, ncol = ncol(Y), nrow = nrow(Y))
    X <- sign(Y) * pmax(abs(Y) - c, myzero, na.rm = TRUE)
    X
    } 

# Test
norm(prox_l1(dat_11, 5), "F")

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

# Test
prox_nuclear(dat_11, 5)
norm(prox_nuclear(x, 5)[[1]], type = "F")

############################################################
############################################################

# is same function for convergence criteria
# is the difference among matrices > noise threshold?
## if TRUE, keep iterating, if FALSE, end

# Compares L1, L2, L3 OR S1, S2
# Need ... for function to handle different number of inputs
# length(varargin) gives the number of function input arguments given in the call
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

# Test
is_same(0.5, x, (x-0.5), (x-0.4))
is_same(0.5, x, (x-0.4), (x-0.7))
is_same(0.5, x, (x-0.4), (x+0.2))

############################################################
############################################################

# loss_lod function (only used in the loss function)

loss_lod <- function(X, D, Delta) {
  # % D is the original data
  # % X is the new thing (L + S)
  # Delta is the LOD
  
  # % Pointwise boolean operation tricks for element-wise updating
  X_lod <- (X - D)     * (D >= 0) +
  # % D>=0 will spit out 0/1 (no/yes)
  # % If D_ij >= 0, then X_lod = (X - D)_ij, else zero
  # Normal loss for >LOD measures (distance from original value)
    
           t(t(X) - Delta) * (D < 0 & t(t(X) > Delta)) +
  # % If D_ij < 0 AND X_ij > Delta, then X_lod = X_ij - Delta, else zero
  # % D is < 0 when < LOD
  # This should be penalized more because D <LOD but (L+S) >LOD (distance from LOD)
    
            X          * (D < 0 & X < 0)
  # % If D_ij < 0 AND X_ij < 0, then X_lod = X, else zero
  
  l <- sum(X_lod^2) / 2
  # % L2 norm
  
  # % Any D_ij < 0 AND X_ij < Delta AND > 0 are treated as equal
  
  # % Minimize discrepancy for valid data
  # % Want to shrink negative things
  l
}

# Test
loss_lod(dat_11, (dat_11-100), 1:10)
loss_lod(dat_11, (dat_11-100), 0)
loss_lod(x, (x-5), delta50)
loss_lod(x, (x-5), 7)
# if delta is scalar, this is the same
# if delta is vector, NOT THE SAME

# DIFFERENT
dat > (1:10)

# SAME
t(t(dat) > 1:10)
t(t(dat) > c(20, 5, 22:29))

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
  
  # Max iteration
  MAX_ITER <- 100
  
  # Convergence Thresholds
  LOSS_THRESH <- 1e-5
  SAME_THRESH <- 1e-4
  
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
    
    # These are all derivatives
    L2_opt1 <- (mu*rho*D     + (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)
    L2_opt2 <- L1 + Z1/rho
    L2_opt3 <- t((mu*rho*Delta + t(((mu + rho)*Z1) - (mu*Z3) + ((mu + rho)*rho*L1) - (mu*rho*S1)))) / ((2*mu*rho) + (rho^2))
    L2_opt4 <- (               (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)
    
    L2_new <- (L2_opt1 * (D >= 0)) +
      # If D >= LOD, use opt1 (Good)
              (L2_opt2 * ((D < 0) & ((L2 + S2) >= 0) & t(t(L2 + S2) <= Delta))) +
      # If D < LOD and new is between 0 and LOD, use opt2 (Good)
              (L2_opt3 * ((D < 0) & t(t(L2 + S2) > Delta))) +
      # If D < LOD and new > LOD use opt3 (Bad)
              (L2_opt4 * ((D < 0) & ((L2 + S2) < 0)))
      # If D < LOD and new < LOD, use opt4 (Bad)
      # % L2_new becomes whichever of the 4 meets the conditions
    
    S2_opt1 <- (mu*rho*D     + (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)
    S2_opt2 <- S1 + (Z3/rho)
    S2_opt3 <- t(((mu*rho*Delta) + t(((mu + rho)*Z3) - (mu*Z1) + ((mu + rho)*rho*S1) - (mu*rho*L1)))) / ((2*mu*rho) + (rho^2))
    S2_opt4 <- (               (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)
    
    S2 <- (S2_opt1 * (D >= 0)) +
          (S2_opt2 * ((D < 0) & ((L2 + S2) >= 0) & t(t(L2 + S2) <= Delta))) +
          (S2_opt3 * ((D < 0) & t(t(L2 + S2) > Delta))) +
          (S2_opt4 * ((D < 0) & ((L2 + S2) < 0)))
    # % For data >LOD, use opt 1
    # % S2 becomes whichever of the 4 meets the conditions
    
    L2 <- L2_new
    # % The code block above takes LOD into account.
    # % The code block commented out below does not take LOD into account
    # %     L2 = (mu*rho*D + (mu+rho)*Z1 - mu*Z3 + (mu+rho)*rho*L1 - mu*rho*S1) / (2*mu*rho+rho^2);
    # %     S2 = (mu*rho*D + (mu+rho)*Z3 - mu*Z1 + (mu+rho)*rho*S1 - mu*rho*L1) / (2*mu*rho+rho^2);
    
    L3 <- pmax(L1 + Z2/rho, 0, na.rm = TRUE)
    # % Non-Negativity constraint!
    
    Z1 <- Z1 + rho*(L1 - L2)
    Z2 <- Z2 + rho*(L1 - L3)
    Z3 <- Z3 + rho*(S1 - S2)
    # % Z accumulate differnces between L and L and between S and S
    
    loss[i] <- nuclearL1 + 
      (lambda*sum(abs(S1))) +
      (mu*loss_lod((L2 + S2), D, Delta)) +
      sum(Z1*(L1 - L2)) +
      sum(Z2*(L1 - L3)) +
      sum(Z3*(S1 - S2)) +
      (rho/2 * (sum((L1-L2)^2) + sum((L1 - L3)^2) + sum((S1 - S2)^2)))
    # % The code block above takes LOD into account.
    
    # % The code block commented out below does not take LOD into account
    # %     loss(i) = nuclearL1 + lambda*sum(sum(abs(S1))) + mu/2*sum(sum((L2+S2-D).^2)) ...
    # %         + sum(sum(Z1.*(L1-L2))) + sum(sum(Z2.*(L1-L3))) + sum(sum(Z3.*(S1-S2))) ...
    # %         + rho/2 * ( sum(sum((L1-L2).^2)) + sum(sum((L1-L3).^2)) + sum(sum((S1-S2).^2)) );

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

# Test
# R results 
m <- nrow(mixture_data)
r_output <- pcp_lod(x, 4/sqrt(m), 10, 0)
lowrank_r <- r_output$L
sparse_r <- r_output$S

# MATLAB results
lowrank_m <- readMat("./Below_LOD/MATLAB/LOD_demo_output/lowrank_lod0.mat") %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  as.matrix()
sparse_m <- readMat("./Below_LOD/MATLAB/LOD_demo_output/sparse_lod0.mat") %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  as.matrix()

# Pretty close!
norm(lowrank_r - lowrank_m, type = "F")
norm(sparse_r - sparse_m, type = "F")

head(lowrank_r)[1:5]
head(lowrank_m)[1:5]

######################################
# 50% <LOD
# MATLAB results
lowrank_m50 <- readMat("./Below_LOD/MATLAB/LOD_demo_output/lowrank_lod50.mat") %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  as.matrix()
sparse_m50 <- readMat("./Below_LOD/MATLAB/LOD_demo_output/sparse_lod50.mat") %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  as.matrix()

# R results 
r_output50 <- pcp_lod(mix_data_lod_50, 4/sqrt(m), 10, delta50)
lowrank_r50 <- r_output50$L
sparse_r50 <- r_output50$S

# same-ish
norm(lowrank_r50 - lowrank_m50, type = "F")
norm(sparse_r50 - sparse_m50, type = "F")

head(lowrank_r50)[1:5]
head(lowrank_m50)[1:5]
