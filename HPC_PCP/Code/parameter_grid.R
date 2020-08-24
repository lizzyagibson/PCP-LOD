
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(pcpr, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(pracma, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(Matrix, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

############################################
# 0. Grab HPC job id
############################################

## read job number from system environment
## This only works if run on cluster!!
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num 

############################################
# 1. Set up parameter grid
############################################

# matrix dim
n <- 2500 # num. rows
m <- 20 # num. cols
r <- 6 # rank

## Vary $\mu$ and $\lambda$
mu_values <- seq(0.01, sqrt(m/(2*log(n*m))), by = 0.05)
lam_values <- seq(1/sqrt(n), 1/sqrt(m), 0.01)

seeds = 1:100

param_grid <- expand_grid(seeds, lam_values, mu_values)

############################################
# 2. Simulate data
############################################

sigma <- 0.3 # std dev for gaussian noise (Z)

# function to simulate a data matrix / chemical mixture...
sim_data <- function(sim_seed, dimr=n, dimc=m, mrank=r, noise_std=sigma) {
  # gaussian noise = Z:
  set.seed(994 + sim_seed)
  Z <- randn(dimr, dimc) * noise_std
  
  # low rank matrix = L:
  set.seed(1996 + sim_seed)
  U <- rand(dimr, mrank)
  set.seed(2998 + sim_seed)
  V <- rand(mrank, dimc)
  L <- U%*%V # ground truth low rank matrix. 
  
  # intermediate tstep o help define sparse matrix S below:
  D <- L + Z
  set.seed(sim_seed)
  
  # sparse matrix = S:
  S <- -D * ((D < 0) * 1) + (rand(n,m)<0.03) * rand(n,m)*1 # Jingkai's method of simulating sparse noise (1st term ensures nonnegativity, second term adds some sparse noise to random entries)
  #S <- (-D + rand(dimr, dimc)) * ((D < 0) * 1) # Lawrence's method of simulating sparse noise (use sparse to ensure matrix is non-neg + add extra noise on top of those same entries)
  #S <- rsparsematrix(nrow = dimr, ncol = dimc, density = .07, rand.x = rand) # for when there is no gaussian (Z) noise added. buggy at the moment throwing fatal error.
  
  # final simulated data matrix = M:
  X <- S + D
  
  list(X = X, L = L, S = S)
}

############################################
# 3. PCP function
############################################

pcp_func <- function(X, lam, mu){ 
  mixture <- pcp_lod(X, lam, mu, 0)
  L <- mixture$L
  S <- mixture$S
  sv_diag <- svd(L)$d # singular values on new low rank matrix
  rank_L <- rankMatrix(L)

  list(L_out = L, SV_L = sv_diag, S_out = S, rank_L = rank_L)
}

############################################
# 4. Loop over grid
############################################

# Loop over $\lambda$, $\mu$ pairs.
range_out <- param_grid[1, ] %>% 
  mutate(sim_out = map(seeds, sim_data)) %>% 
  unnest_wider(sim_out) %>% 
  mutate(pcp_out = pmap(list(X, lam_values, mu_values),
                        pcp_func)) %>% 
  unnest_wider(pcp_out)

############################################
# 5. Evaluate Sparsity
############################################

cells <- n*m

prop_dense <- function(sparse){
  as_tibble(sparse) %>% 
    pivot_longer(V1:V20) %>% 
    mutate(binary = ifelse(value == 0, 0, 1)) %>% 
    summarize(density = sum(binary)/cells) %>% 
    pull(density)
  }

sparse_out <- range_out %>% 
  mutate(density_S = map(S_out, prop_dense)) %>% 
  unnest(density_S)

############################################
# 6. Evaluate Error
############################################

get_L_error <- function(L, L_out) {
  norm((L - L_out), type = "F")/norm(L, type = "F")
  }

all_out <- sparse_out %>% 
  mutate(L_rel_error = map2(L, L_out, get_L_error)) %>% 
  unnest(L_rel_error)

save(all_out, file = paste0("param_out_", job_num, ".RDA"))







