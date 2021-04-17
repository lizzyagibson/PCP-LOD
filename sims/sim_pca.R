# PCA on simulated datasets

library(tidyverse)
library(GGally)
library(PCPhelpers)
library(pcpr)
library(factoextra)

load("./sims/sim_lod.RDA")
sim_lod

get_pca <- function (sim) {
  # Run PCA centered, not scaled
  pca_out <- prcomp(sim)
  #loadings
  rot <- pca_out$rotation
  # scores
  ex <- pca_out$x
  # singular values
  sv <- pca_out$sdev

  rank = 4
  # Cut scores and patterns to rank
  rotations <- as_tibble(rot[, 1:rank])
  scores <- if (rank == 1) {matrix(ex[, 1:rank], nrow = nrow(sim))} else {ex[, 1:rank]}
  # Predicted values
  pred <- scores %*% t(rotations) + kronecker(matrix(1, 500, 1), t(apply(sim, 2, mean)))

  # Explain >=80% of var
  pve <- sv^2/sum(sv^2)
  rank <- 0
  for (i in 1:length(sv)) {
    if (sum(pve[1:i]) >= 0.8) {
      rank <- i
      break
    }}
  
  return(list(rotations = rotations, scores = scores, pred = pred, rank = rank))
}

# rotations are right singular vectors

# Run PCA ####
sim_pca = sim_lod %>% 
          mutate(lod_sqrt2_mat = map(lod_sqrt2_mat, function(x) apply(x, 2, function(y) y/sd(y))),
                 pca_out      = map(lod_sqrt2_mat, get_pca),
                 pca_loadings = map(pca_out, function(x) x[[1]]),
                 pca_scores   = map(pca_out, function(x) x[[2]]),
                 pca_pred     = map(pca_out, function(x) x[[3]]),
                 pca_rank     = map(pca_out, function(x) x[[4]]),
                 pca_error    = map2(chem, pca_pred, function(x,y)
                                                      norm(x-y,"F")/norm(x,"F"))) %>% 
  unnest(c(pca_rank, pca_error))

# Quick summary ####
summary(sim_pca$pca_error)
summary(sim_pca$pca_rank)
sum(sim_pca$pca_rank == 4)/600

fviz_eig(prcomp(sim_lod$lod_sqrt2_mat[[3]])) 

# Get metrics ####
pca_metrics = sim_pca %>% 
  mutate(method = "PCA") %>% 
  arrange(seed, chem) %>% 
  mutate(all_relerr   = map2(chem, pca_pred, function(x,y) norm(x-y,"F")/norm(x,"F"))) %>% 
  unnest(c(all_relerr)) %>% 
  mutate(mask = map(lod_neg1_mat, function(x) x != -1),
         above_relerr = pmap(list(mask, chem, pca_pred), 
                               function(mask, x,y) norm(mask*x-mask*y,"F")/norm(mask*x,"F")),
         below_relerr = pmap(list(mask, chem, pca_pred), 
                               function(mask, x,y) norm((!mask)*x-(!mask)*y,"F")/norm((!mask)*x,"F"))) %>% 
  unnest(c(above_relerr, below_relerr)) 

all(pca_metrics$all_relerr == pca_metrics$pca_error)

sim_pca %>% 
  group_by(chemicals, lim) %>% 
  summarize(qs = quantile(pca_error), prop = seq(0, 1, 0.25)) %>% 
  pivot_wider(names_from = prop,
              values_from = qs)

pca_combine = pca_metrics %>% 
  dplyr::select(seed, chemicals, lim, method, all_relerr, above_relerr, below_relerr) %>% 
  mutate(LS = "PCA")

