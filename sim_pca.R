# PCA on simulated datasets

library(tidyverse)
library(GGally)
library(PCPhelpers)
library(pcpr)

load("./sims/sim_lod.RDA")

get_pca <- function (sim) {
  # Run PCA centered, not scaled
  pca_out <- prcomp(sim)
  #loadings
  rot <- pca_out$rotation
  # scores
  ex <- pca_out$x
  # singular values
  sv <- pca_out$sdev
  
  # Explain >=80% of var
  pve <- sv^2/sum(sv^2)
  rank <- 0
  for (i in 1:length(sv)) {
    if (sum(pve[1:i]) >= 0.8) {
      rank <- i
      break
    }}
  
  # Cut scores and patterns to rank
  rotations <- as_tibble(rot[, 1:rank])
  scores <- if (rank == 1) {matrix(ex[, 1:rank], nrow = nrow(sim))} else {ex[, 1:rank]}
  # Predicted values
  pred <- scores %*% t(rotations) + kronecker(matrix(1, 500, 1), t(apply(sim, 2, mean)))

  
  return(list(rotations = rotations, scores = scores, pred = pred, rank = rank))
}

sim_pca = sim_lod %>% 
          mutate(lod_sqrt2_mat = map(lod_sqrt2_mat, function(x) apply(x, 2, function(y) y/sd(y))),
                 pca_out      = map(lod_sqrt2_mat, get_pca),
                 pca_loadings = map(pca_out, function(x) x[[1]]),
                 pca_scores   = map(pca_out, function(x) x[[2]]),
                 pca_pred     = map(pca_out, function(x) x[[3]]),
                 pca_rank     = map(pca_out, function(x) x[[4]]),
                 pca_error    = map2(chem, pca_pred, function(x,y)
                                                      norm(x-y,"F")/norm(x,"F")))

summary(unlist(sim_pca$pca_error))

sum(unlist(sim_pca$pca_rank) == 4)/600

sim_pca$pca_loadings[[1]] %>% 
  mutate(id = fct_inorder(str_c("Chem ", 1:nrow(.)))) %>% 
  pivot_longer(PC1:PC8) %>% 
  ggplot(aes(x = id, y = name)) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = "RdYlBu") + theme_test() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

sim_pca$pca_loadings[[1]] %>% 
  mutate(id = fct_inorder(str_c("Chem ", 1:nrow(.)))) %>% 
  pivot_longer(PC1:PC8) %>% 
  ggplot(aes(x = id, y = value)) +
  geom_segment( aes(x=id, xend=id, y=0, yend=value), color="grey") +
  geom_point( color="orange", size=1) +
  facet_wrap(.~name) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1))

