# PCA on simulated datasets

# Load functions
source("Sims/functions.R")

# Load data
load("./sims/sim_lod.RDA")
sim_lod

# Run this function on all simulations
get_pca <- function (sim) {
  # Run PCA 
  pca_out <- prcomp(sim)
  
  #loadings
  rot <- pca_out$rotation
  # scores
  ex <- pca_out$x
  # singular values
  sv <- pca_out$sdev

  # Cut scores and patterns to this rank
  # Explain >=80% of var
  pve <- sv^2/sum(sv^2)
  rank <- 0
  for (i in 1:length(sv)) {
    if (sum(pve[1:i]) >= 0.8) {
      rank <- i
      break
    }}
  
  rotations <- as_tibble(rot[, 1:rank])
  scores <- if (rank == 1) {matrix(ex[, 1:rank], nrow = nrow(sim))} else {ex[, 1:rank]}
  # Predicted values
  pred <- scores %*% t(rotations) + kronecker(matrix(1, 500, 1), t(apply(sim, 2, mean)))
  # this adds back the mean
  
  return(list(rotations = rotations, scores = scores, pred = pred, rank = rank))
}

# Run PCA ####
sim_pca = sim_lod %>% 
          mutate(pca_out      = map(lod_sqrt2_mat, get_pca),
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
sum(sim_pca$pca_rank == 4)/1800

# Get metrics ####
pca_metrics = pcp_out_re %>% 
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

