# CV

library(tidyverse)
library(GGally)
library(PCPhelpers)
library(pcpr)
library(patchwork)

# 4 patterns
# 2 mixture sizes, 16 & 48
# 1 sample size = 500
# 3 proportions <LOD, 25, 50, 75

load("./sims/sim_lod.RDA")
sim_lod

# Run PCP-LOD on all sims!
sim_pcp_out = sim_lod %>% 
              filter(lim != 0.75) %>% 
              mutate(pcp_out = map2(lod_neg1_mat, lod, function(x,y)
                        root_pcp_noncvx_nonnegL_na_lod(D = x, MAX_ITER = 20000,
                                                  lambda = 1/sqrt(nrow(x)), mu = sqrt(ncol(x)/2), 
                                                  r = 4, LOD = y,
                                                          verbose = TRUE)))
#save(sim_pcp_out, file = "./sims/sim_pcp_25_50.rda")
sim_lod_75 = sim_lod %>% 
  filter(lim == 0.75)
#save(sim_lod_75, file = "./sims/sim_lod_75.rda")

# sim_pcp_out_75 = sim_lod_75 %>% 
#                   mutate(pcp_out = map2(lod_neg1_mat, lod, function(x,y)
#                          root_pcp_noncvx_nonnegL_na_lod(D = x, MAX_ITER = 100000,
#                                    lambda = 1/sqrt(nrow(x)), mu = sqrt(ncol(x)/2), 
#                                    r = 4, LOD = y,
#                                    verbose = TRUE)))
# pcp_out_75 = tibble()
# for (i in 1:200) {
#   load(paste0("/ifs/scratch/msph/ehs/eag2186/pcp/lod/sim_pcp_out_75_", i, ".rda"))
#   pcp_out_75 <- bind_rows(pcp_out_75, sim_pcp_out_75)
#        print(i)
# }
load("./sims/pcp_out_75.rda")

pcp_out = sim_pcp_out %>% 
            bind_rows(., pcp_out_75) %>% 
        arrange(seed, chem) %>% 
            mutate(L = map(pcp_out, function(x) x$L),
                   S = map(pcp_out, function(x) x$S),
                   LS = map2(L,S,`+`)) %>% 
  mutate(mse_l = map2(sim, L, function(x,y) mean((x-y)^2)),
         mse   = map2(sim, LS, function(x,y) mean((x-y)^2)),
         relerr_l = map2(sim, L, function(x,y) norm(x-y,"F")/norm(x,"F")),
         relerr   = map2(sim, LS, function(x,y) norm(x-y,"F")/norm(x,"F")))

pcp_out %>% unnest(c(mse_l, mse, relerr_l, relerr))

summary(unlist(pcp_out$mse))
summary(unlist(pcp_out$mse_l))

summary(unlist(pcp_out$relerr))
summary(unlist(pcp_out$relerr_l))

all(sim_lod$sim[[600]] == pcp_out$sim[[600]])
