# PCA on simulated datasets

# Load functions
source("./functions.R")

# Load data
load("./Sims/Sim Data/sim_lod.RDA")
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

# SVD ####
# need to rerrange eigen vectors to best align with those of measured data
pca_svd_out = sim_pca %>% 
  mutate(svd_chem = map(chem, svd),
         svd_pca  = map(lod_sqrt2_mat, svd),
         svd_chem_left  = map(svd_chem, function(x) x$u),
         svd_chem_right = map(svd_chem, function(x) x$v),
         svd_pca_left  = map(svd_pca, function(x) x$u),
         svd_pca_right = map(svd_pca, function(x) x$v))

# Rearrange vectors to be closest to original SVD results
pca_svd_re = pca_svd_out %>% 
             mutate(svd_pca_left  = map2(svd_chem_left,  svd_pca_left,  factor_correspondence),
                    svd_pca_right = map2(svd_chem_right, svd_pca_right, factor_correspondence))

save(pca_svd_re, file = "./Sims/Sim Data/pca_svd_re.rda")  
# load("./Sims/Sim Data/pca_svd_re.rda")  

# Get metrics ####
pca_svd_metrics = pca_svd_re %>% 
  mutate(method = "PCA") %>% 
  mutate(left   = map2(svd_chem_left, svd_pca_left, function(x,y) norm(x-y,"F")/norm(x,"F")),
         right   = map2(svd_chem_right, svd_pca_right, function(x,y) norm(x-y,"F")/norm(x,"F"))) %>% 
  select(-c(true_patterns,  true_scores,  chem,   sim, pca_loadings,  pca_scores,   pca_pred,
            lod,   lod_neg1_mat, lod_sqrt2_mat, svd_chem, svd_pca, sparsity, pca_out,
            svd_chem_left, svd_pca_left, svd_pca_right, svd_chem_right)) %>% 
  unnest(c(left, right))

pca_metrics = sim_pca %>% 
  mutate(method = "PCA") %>% 
  arrange(seed, chem) %>% 
  mutate(all_relerr   = map2(chem, pca_pred, function(x,y) norm(x-y,"F")/norm(x,"F"))) %>% 
  unnest(c(all_relerr)) %>% 
  mutate(mask = map(lod_neg1_mat, function(x) x != -1),
         above = pmap(list(mask, chem, pca_pred), 
                               function(mask, x,y) norm(mask*x-mask*y,"F")/norm(mask*x,"F")),
         below = pmap(list(mask, chem, pca_pred), 
                               function(mask, x,y) norm((!mask)*x-(!mask)*y,"F")/norm((!mask)*x,"F"))) %>% 
  select(-c(true_patterns, true_scores,  chem,   sim, sparsity,
            pca_out, pca_loadings, pca_scores, pca_pred, pca_error, pca_rank,
            lod, lod_neg1_mat, lod_sqrt2_mat, mask)) %>% 
  unnest(c(above, below))

# Save metrics
save(pca_metrics, file = "./Sims/Sim Data/pca_metrics.rda")
save(pca_svd_metrics, file = "./Sims/Sim Data/pca_svd_metrics.rda")
