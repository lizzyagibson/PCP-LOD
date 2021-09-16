# CV
# Load packages
source("./functions.R")

# 4 patterns
# 2 mixture sizes, 16 & 48
# 1 sample size = 500
# 3 proportions <LOD, 25, 50, 75

# Load data
load("./Sims/Sim Data/sim_lod.RDA")
sim_lod

# Run PCP-LOD ####
# Run PCP-LOD on all sims!
sim_pcp_out = sim_lod %>% 
              mutate(pcp_out = map2(lod_neg1_mat, lod, function(x,y)
                        root_pcp_noncvx_nonnegL_na_lod(D = x, MAX_ITER = 20000,
                                                  lambda = 1/sqrt(nrow(x)), mu = sqrt(ncol(x)/2),
                                                  r = 4, LOD = y)), # r = 4 because that's what cross validation chose
                      L = map(pcp_out, function(x) x$L))
# Jaime: running lines 16 to 21 I get the following error, Error: Problem with `mutate() x object 'L' not found
# Lizzy: we fixed this with the updated pcpr package

save(sim_pcp_out, file = "./Sims/Sim Data/sim_pcp_out.rda")
# load("./Sims/Sim Data/sim_pcp_out.rda")
sim_pcp_out

# SVD ####
pcp_svd_out = sim_pcp_out %>% 
  mutate(svd_chem = map(chem, svd),
         svd_pcp  = map(L, svd),
         svd_chem_left  = map(svd_chem, function(x) x$u),
         svd_chem_right = map(svd_chem, function(x) x$v),
         svd_pcp_left  = map(svd_pcp, function(x) x$u),
         svd_pcp_right = map(svd_pcp, function(x) x$v))

# Rearrange vectors to be closest to original SVD results
# This step takes a long time, I let it run overnight, you can just load the rda file below
pcp_svd_re = pcp_svd_out %>% 
  mutate(svd_pcp_left  = map2(svd_chem_left,  svd_pcp_left,  factor_correspondence),
         svd_pcp_right = map2(svd_chem_right, svd_pcp_right, factor_correspondence))

save(pcp_svd_re, file = "./Sims/Sim Data/pcp_svd_re.rda")  
# load("./Sims/Sim Data/pcp_svd_re.rda")  

# Get metrics ####
pcp_svd_metrics = pcp_svd_re %>% 
  mutate(method = "PCP-LOD") %>% 
  mutate(left   = map2(svd_chem_left, svd_pcp_left, function(x,y) norm(x-y,"F")/norm(x,"F")),
         right  = map2(svd_chem_right, svd_pcp_right, function(x,y) norm(x-y,"F")/norm(x,"F"))) %>% 
  select(-c(true_patterns, true_scores,  chem,   sim,
            lod, lod_neg1_mat, lod_sqrt2_mat, svd_chem, svd_pcp, pcp_out,
            L, svd_chem_left, svd_pcp_left, svd_pcp_right, svd_chem_right)) %>% 
  unnest(c(left, right))

# Get metrics ####
pcp_metrics = sim_pcp_out %>% 
  mutate(method = "PCP-LOD") %>% 
  arrange(seed, chem) %>% 
  mutate(all_relerr   = map2(chem, L, function(x,y) norm(x-y,"F")/norm(x,"F"))) %>% 
  unnest(c(all_relerr)) %>% 
  mutate(mask = map(lod_neg1_mat, function(x) x != -1),
         above = pmap(list(mask, chem, L), 
                               function(mask, x,y) norm(mask*x-mask*y,"F")/norm(mask*x,"F")),
         below = pmap(list(mask, chem, L), 
                               function(mask, x,y) norm((!mask)*x-(!mask)*y,"F")/norm((!mask)*x,"F"))) %>% 
  select(-c(true_patterns, true_scores,  chem, sim, mask,
            lod, lod_neg1_mat, lod_sqrt2_mat, pcp_out, L)) %>% 
  unnest(c(all_relerr, above, below)) 

# Save metrics
save(pcp_metrics, file = "./Sims/Sim Data/pcp_metrics.rda")
save(pcp_svd_metrics, file = "./Sims/Sim Data/pcp_svd_metrics.rda")


