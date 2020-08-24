
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(pcpr, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(pracma, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(Matrix, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

############################################
# 0. Grab HPC job id
############################################

## read job number from system environment
## This only works if run on cluster!!
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))

############################################
# 1. Read Boston data
############################################

## Real Boston data

# Read air pollution data
mixture <- readMat("Data/mixtures_data.mat")

mix <- as.data.frame(mixture) %>% as_tibble() %>% 
  select(Al, As, Ba, bc, Br, Ca, Cl,
         Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
         Ti,  V, Zn) %>% 
  drop_na(.) %>% 
  as.matrix(.)

nn <- nrow(mix)
mm <- ncol(mix)

#This is the default
lam_mix = 1/sqrt(nn)
mu_mix = sqrt(mm/(2*log(nn*mm)))

############################################
# 2. Set up parameter grid
############################################

# matrix dim
n <- 2500 # num. rows
m <- 20 # num. cols
r <- 6 # rank

## Vary $\mu$ and $\lambda$
mu_values <- seq(0.01, mu_mix, by = 0.05)
lam_values <- seq(1/sqrt(n), 1/sqrt(m), 0.01)

param_grid <- expand_grid(lam_values, mu_values)

############################################
# 3. PCP function
############################################

pcp_func <- function(lam, mu){ 
  mixture <- pcp_lod(mix, lam, mu, 0)
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
range_out <- param_grid[job_num, ] %>% 
  mutate(pcp_out = pmap(list(lam_values, mu_values),
                        pcp_func)) %>% 
  unnest_wider(pcp_out)

############################################
# 5. Evaluate Sparsity
############################################

cells <- nn*mm

prop_dense <- function(sparse){
  as_tibble(sparse) %>% 
    pivot_longer(Al:Zn) %>% 
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

get_total_error <- function(L_out, S_out) {
  norm((mix - L_out - S_out), type = "F")/norm(mix, type = "F")
}

all_out <- sparse_out %>% 
  mutate(total_rel_error = map2(L_out, S_out, get_total_error)) %>% 
  unnest(total_rel_error)

all_out

save(all_out, file = paste0("boston_param_out_", job_num, ".RDA"))












