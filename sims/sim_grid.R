# Create simulated datasets

library(tidyverse)
library(LearnBayes) # this has Dirchlet dist
library(GGally)
library(PCPhelpers)
library(pcpr)
library(pracma)
library(grid)
library(factoextra)
library(RColorBrewer)

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

p = as_tibble(patterns)
colnames(p) = 1:ncol(p)
# pdf("./sims/loadings_plot.pdf", width = 15, height = 15)
p %>%
  mutate(Pattern = 1:nrow(.),
         Pattern = str_c("Pattern ", Pattern)) %>%
  pivot_longer(1:16) %>%
  mutate(name = fct_inorder(str_remove(name, "V"))) %>%
  ggplot(aes(x = name, y = value)) +
  geom_col(fill = "orange", color = "black") +
  facet_wrap(.~Pattern) +
  labs(y = "Simulated Loadings", x = "Simulated Chemicals") +
  theme_light(base_size = 45) +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 20),
        #axis.ticks.x = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size = 45, colour = 'black'))
dev.off()

brewer.pal(8, "Spectral")
display.brewer.pal(8, "Spectral")

p %>%
  mutate(Pattern = 1:nrow(.),
         Pattern = str_c("Pattern ", Pattern)) %>%
  pivot_longer(1:16) %>%
  mutate(name = fct_inorder(str_remove(name, "V"))) %>%
  ggplot(aes(x = name, y = value, fill = Pattern)) +
  geom_col() +
  labs(y = "Simulated Loadings", x = "Simulated Chemicals", fill = "") +
  scale_fill_manual(values = c("#F46D43", "#66C2A5", "#FDAE61", "#2B83BA")) +
  theme_light(base_size = 20) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position = "bottom",
        #axis.ticks.x = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size = 45, colour = 'black'))

# 100 random samples from each data generating process
seed = 1:100
chemicals = c(16, 48)

# get every combination
pattern_comb = expand_grid(seed, chemicals)

# simulate patterns
patterns_iter <- pattern_comb %>% 
  mutate(true_patterns = map2(seed, chemicals, create_4patterns))

# Simulate Scores
create_scores <- function (seed) {
  n = 500
  r = 4
  
  seed = 19888891 + seed
  set.seed(seed)
  # independent draws from standard log-normal
  scores <- matrix(exp(rnorm(n*r, mean = 1)), nrow = n, ncol = r)
  return(as.matrix(scores))
}

pdf("./sims/figures/scores_plot.pdf")
scores %>%
  as_tibble() %>% 
  mutate(Pattern = 1:nrow(.),
         Pattern = str_c("Pattern ", Pattern)) %>%
  pivot_longer(1:4) %>%
  mutate(name = fct_inorder(str_remove(name, "V"))) %>%
  ggplot(aes(x = value, y = ..density..)) +
  geom_histogram(fill = "orange", color = "black") +
  labs(y = "Density", x = "Simulated Scores") +
  theme_light(base_size = 45) +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(size = 30),
        axis.text.x = element_text(size = 20))
dev.off()

scores_iter <- patterns_iter %>% 
  mutate(true_scores = map(seed, create_scores))

# Simulate Chemical Exposures
# matrix multiply scores times loadings
sim_iter <- scores_iter %>% 
  mutate(chem = map2(true_scores, true_patterns, `%*%`))

# Simulate Noise
add_noise <- function (seed, chem, sd) {
  n = nrow(chem)
  p = ncol(chem)
  noise <- matrix(NA, nrow = n, ncol = p)
  
  seed = 19888891 + seed
  set.seed(seed)

  # add noise from normal dist, mean = 0
  for (i in 1:p) {
    noise[,i] <- (rnorm(n, sd = sd)) # was 4 before
  }
  
  # if negative, push to zero
  sim = pmax(chem + noise, 0)
  colnames(sim) = str_c("chem_", str_pad(1:ncol(sim), 2, pad = "0"))
  sim
}

# add noise
sim_iter <- sim_iter %>% 
  mutate(sim_5 = pmap(list(seed, chem), function(x,y) add_noise(x,y,5)),
         sim_1 = pmap(list(seed, chem), function(x,y) add_noise(x,y,1)))

sum(sim_iter$sim_5[[1]] == 0)/(500*16)

# sort(1/apply(sim_iter$chem[[1]], 2, sd))
# heatmaply::heatmaply(cor(sim_iter$sim[[1]]))

add_sparse <- function (seed, sim) {
  n = nrow(sim)
  p = ncol(sim)

  seed = 19888891 + seed
  set.seed(seed)

  S <- (rand(n,p)<0.05) * (rand(n,p)*(max(sim)/4))
  # adds some sparse noise to random entries
  # approx 5% sparse events

  simS = sim + S

  list(simS = simS, S = S)
}

# apply(sim,2,max)
# apply(simS,2,max)

# add sparse
sim_iter <- sim_iter %>%
 mutate(sparse_out = map2(seed, sim_1, add_sparse),
        sim_sparse = map(sparse_out, function(x) x$simS),
        sparsity = map(sparse_out, function(x) x$S))

# Next, we subject simulated chemicals to an LOD:
# Using corrupt_mat from pcphelpers
lim_range = expand_grid(pattern_comb, lim = c(0.25,0.5,0.75))

sim_lod = sim_iter %>% 
        dplyr::select(-sparse_out) %>% 
        pivot_longer(sim_5:sim_sparse,
                     values_to = "sim") %>% 
        left_join(., lim_range) %>% 
          mutate(lod = map2(sim, lim, function(x,y)
                            as.vector(apply(x, 2, quantile, probs = y))),
                 lod_neg1_mat = map2(sim, lim, function(x,y) 
                                 corrupt_mat(x, cols = 1:ncol(x), limit=y, fill="-1")),
                 lod_sqrt2_mat = map2(sim, lim, function(x,y) 
                               corrupt_mat(x, cols = 1:ncol(x), limit=y, fill="sqrt2")))

# save nested dataframe 
# save(sim_lod, file = "./sims/sim_lod.RDA")

get_lower_tri <-function(x){
  x[upper.tri(x)] <- NA
  return(x)
}

cormat <- as.data.frame(get_lower_tri(cor(sim_lod$sim[[1]])[c(5,6,3,7:10,12,11,13:15,4,16,1,2),
                                                                    c(5,6,3,7:10,12,11,13:15,4,16,1,2)])) %>% 
  rownames_to_column(var = "Chem") %>% 
  as_tibble() %>% 
  mutate(Chem = fct_inorder(Chem),
         Chem = as_factor(1:nrow(.))) %>% 
  pivot_longer(chem_05:chem_02) %>% 
  mutate(name = fct_inorder(name),
         name = as_factor(rep(1:16,16)))

corplot = cormat %>%
  mutate(outline = ifelse(is.na(value), FALSE, TRUE))
corplot$outline[!corplot$outline] <- NA

#pdf("./sims/sim_corr.pdf")
corplot %>% 
  ggplot(aes(x = Chem, y = name, fill = value)) + 
  geom_tile() +
  geom_tile(data = corplot[!is.na(corplot$outline), ], aes(color = outline), size = 0.75) +
  scale_color_manual(guide = FALSE, values = c(`TRUE` = "black")) +
  labs(x = "Simulated Chemicals", y = "", fill = "") +
  theme_light(base_size = 25) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.075,0.85),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        axis.text.y = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size = 25, colour = 'black')) +
  scale_fill_distiller(palette="Spectral",
                       limits=c(-1,1),
                       na.value = 'white')+
  #scale_color_manual(values = c("black", "white")) +
  guides(color = FALSE)
#dev.off()

  