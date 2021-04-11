# Create simulated datasets

library(tidyverse)
library(LearnBayes) # this has Dirchlet dist
library(GGally)
library(PCPhelpers)
library(pcpr)

# 4 patterns
# 2 mixture sizes, 16 & 48
# 1 sample size = 500
# 3 proportions <LOD, 25, 50, 75

# Simulate Patterns
create_4patterns <- function (seed, chemicals) { 
  seed = 19888891 + seed
  set.seed(seed)
  
  group <- chemicals / 4 # because 4 patterns
  half  <- group/2 # because each group will be half distinct/half overlapping
  
  pat3=c(0,0,1,1) # chemical loads on two patterns and loadings sum to 1
  pat2=c(0,1,1,0)
  pat1=c(1,1,0,0)
  pat4=c(1,0,0,1)
  
  distinct = c() # half of all chemicals load only on 1 pattern
  for (i in 1:half) {
      distinct = cbind(distinct,diag(1,4))
    }
  distinct = distinct[,order(distinct[1,], distinct[2,], distinct[3,], decreasing=T)]
  
  # other half of chemicals load on 2 patterns
  overlap= cbind(t(rdirichlet(half, pat1)), t(rdirichlet(half, pat2)),
                     t(rdirichlet(half, pat3)), t(rdirichlet(half, pat4)))
  
  patterns = cbind(distinct[,1:half], overlap[,1:half], 
                    distinct[,(half+1):(2*half)], overlap[,(half+1):(2*half)],
                    distinct[,(2*half+1):(3*half)], overlap[,(2*half+1):(3*half)],
                    distinct[,(3*half+1):(4*half)], overlap[,(3*half+1):(4*half)])

  return(patterns)
}

# 100 random samples from each data generating process
seed = 1:100
chemicals = c(16, 48)

# get every combination
pattern_comb = expand_grid(seed, chemicals)

# simulate patterns
patterns_iter <- pattern_comb %>% 
  mutate(true_patterns = map2(seed, chemicals, create_4patterns))

patterns_iter$true_patterns[[1]]
# Simulate Scores
create_scores <- function (seed) {
  n = 500
  r = 4
  
  seed = 19888891 + seed
  set.seed(seed)
  # independent draws from standard log-normal
  scores <- matrix(exp(rnorm(n*r)), nrow = n, ncol = r)
  return(as.matrix(scores))
}

scores_iter <- patterns_iter %>% 
  mutate(true_scores = map(seed, create_scores))

# Simulate Chemical Exposures
# matrix multiply scores times loadings
sim_iter <- scores_iter %>% 
  mutate(chem = map2(true_scores, true_patterns, `%*%`))

# Simulate Noise
add_noise <- function (seed, chem) {
  n = nrow(chem)
  p = ncol(chem)
  noise <- matrix(NA, nrow = n, ncol = p)
  
  seed = 19888891 + seed
  set.seed(seed)

  # add noise from normal dist, mean = 0
  for (i in 1:p) {
    noise[,i] <- (rnorm(n, mean = 0, sd = 1))
  }
  
  # if negative, push to zero
  sim = pmax(chem + noise, 0)
  
  colnames(sim) = str_c("chem_", str_pad(1:ncol(sim), 2, pad = "0"))
  
  sim
}

# add noise
sim_iter <- sim_iter %>% 
  mutate(sim = pmap(list(seed, chem), add_noise))

# N(0,1) noise is between 45% and 76% of underlying noise for each chemical
sort(1/apply(sim_iter$chem[[1]], 2, sd))

corr = cor(sim_iter$sim[[1]])
corr[corr == 1] = 0
apply(corr, 2, max)
heatmaply::heatmaply(cor(sim_iter$sim[[1]]))
ggcorr(sim_iter$sim[[4]])

# Next, we subject simulated chemicals to an LOD:
# Using corrupt_mat from pcphelpers
lim_range = expand_grid(pattern_comb, lim = c(0.25,0.5,0.75))

sim_lod = sim_iter %>% 
          left_join(., lim_range) %>% 
          mutate(lod = map2(sim, lim, function(x,y)
                            as.vector(apply(x, 2, quantile, probs = y))),
                 lod_neg1_mat = map2(sim, lim, function(x,y) 
                                 corrupt_mat(x, cols = 1:ncol(x), limit=y, fill="-1")),
                 lod_sqrt2_mat = map2(sim, lim, function(x,y) 
                               corrupt_mat(x, cols = 1:ncol(x), limit=y, fill="sqrt2")))
         
# save nested dataframe 
save(sim_lod, file = "./sims/sim_lod.RDA")



