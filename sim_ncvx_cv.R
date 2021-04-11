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

# CV 1 example for each mixture size/<LOD combo

grid.rank <- tibble(r = 1:10)
sim_cv = sim_lod %>% 
         filter(seed == 1) %>%
         mutate(n = map(sim, nrow),
               p = map(sim, ncol)) %>% 
         unnest(c(n,p)) %>% 
         mutate(mu = sqrt(p/2),
               lambda = 1/sqrt(n),
               noncvx_search = map2(lod_neg1_mat, lod, function(x,y)
                              grid_search_cv(mat = x,
                                             pcp_func = root_pcp_noncvx_nonnegL_na_lod,
                                             grid_df = grid.rank,
                                             lambda = 1/sqrt(500),  mu = sqrt(ncol(x)/2), runs = 5,
                                             cores = 2, LOD = y)))

# save(sim_cv, file = "./sims/sim_cv.RDA")

sim_cv$noncvx_search[[1]]$formatted %>% arrange(desc(value))
sim_cv$noncvx_search[[2]]$formatted %>% arrange(desc(value))
sim_cv$noncvx_search[[3]]$formatted %>% arrange(desc(value))
sim_cv$noncvx_search[[4]]$formatted %>% arrange(desc(value))
sim_cv$noncvx_search[[5]]$formatted %>% arrange(desc(value))
sim_cv$noncvx_search[[6]]$formatted %>% arrange(desc(value))



