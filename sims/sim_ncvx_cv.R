# Cross validate to choose PCP rank

# Get functions
source("Sims/functions.R")

# 4 patterns
# 2 mixture sizes, 16 & 48
# 1 sample size = 500
# 3 proportions <LOD, 25, 50, 75

# Load data
load("./Sims/Sim Data/sim_lod.RDA")
sim_lod

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
                              grid_search_cv(mat = x, # grid_search_cv is in PCPhelpers
                                             pcp_func = root_pcp_noncvx_nonnegL_na_lod,
                                             grid_df = grid.rank,
                                             lambda = 1/sqrt(500),  mu = sqrt(ncol(x)/2), 
                                             runs = 2,
                                             # I ran runs = 100, but that took forever, runs = 2 still takes a while
                                             cores = 2, LOD = y)))

# save(sim_cv, file = "./Sims/Sim Data/sim_cv.RDA")
load("./Sims/Sim Data/sim_cv.RDA")

which.min(sim_cv$noncvx_search[[1]]$formatted$value)
which.min(sim_cv$noncvx_search[[2]]$formatted$value)
which.min(sim_cv$noncvx_search[[3]]$formatted$value)
which.min(sim_cv$noncvx_search[[4]]$formatted$value)
which.min(sim_cv$noncvx_search[[5]]$formatted$value)
which.min(sim_cv$noncvx_search[[6]]$formatted$value)
# CV chose rank 4 for all simulations


