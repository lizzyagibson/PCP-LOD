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

# save(sim_pcp_out, file = "./Sims/Sim Data/sim_pcp_out.rda")
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

# save(pcp_svd_re, file = "./Sims/Sim Data/pcp_svd_re.rda")  
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

##################################################################

# Describe sparsity ####
# sparsity (from sim_lod) has added sparse event

# Across PCP-LOD solutions, between 2% and 10% of S entries were non-sparse. We found decreasing sparsity as the
# proportion < LOD increased, with 3% (IQR: 2%, 4%), 6% (IQR: 4%, 7%), and 7% (IQR: 3%, 8%) unique events, on average, 
# found in simulations with 25%, 50%, and 75% < LOD, respectively
sim_sparse = sim_pcp_out %>% 
  # Don't filter to only sparse events
  # filter(name == "sim_sparse") %>% 
  mutate(S = map(pcp_out, function(x) x$S)) %>% 
  dplyr::select(-c(true_patterns,true_scores,lod,lod_sqrt2_mat, pcp_out)) %>% 
  mutate(lod_mask = map(lod_neg1_mat, function(x) x != -1),
         resid = pmap(list(lod_mask, sim, L), function(lod_mask, sim, L) (sim*lod_mask - L*lod_mask)),
         # sim is truth after noise was added
         # take lod into account
         sparse_thresh = map(resid, function(x) apply(x, 2, sd)*2),
         S_thresh_p = map2(S, sparse_thresh, function(S, sparse_thresh) (S > sparse_thresh)*1),
         S_thresh_n = map2(S, sparse_thresh, function(S, sparse_thresh) (S < -sparse_thresh)*-1),
         S_thresh = map2(S_thresh_p, S_thresh_n, function(x,y) abs(x) + abs(y)))

sparse_overall = sim_sparse %>% 
  mutate(non_sparse = map(S_thresh, function(x) sum(x)/(nrow(x)*ncol(x)))) %>% 
  unnest(non_sparse) %>% 
  group_by(lim) %>% 
  summarize(quant25 = quantile(non_sparse, probs = .25), 
            quant50 = quantile(non_sparse, probs = .5),
            quant75 = quantile(non_sparse, probs = .75)) %>% 
  mutate_all(round, 2) %>% 
  mutate(name = "Overall") %>% select(name, everything())

sparse_size = sim_sparse  %>% 
  mutate(non_sparse = map(S_thresh, function(x) sum(x)/(nrow(x)*ncol(x)))) %>% 
  unnest(non_sparse) %>% 
  group_by(name, lim) %>% 
  summarize(quant25 = quantile(non_sparse, probs = .25), 
            quant50 = quantile(non_sparse, probs = .5),
            quant75 = quantile(non_sparse, probs = .75)) %>% 
  mutate_all(round, 2)

sparse_percent = bind_rows(sparse_size, sparse_overall) %>% 
  pivot_wider(names_from = lim,
              values_from = c(quant25, quant50, quant75)) %>% 
  ungroup() %>% select(c(1, 2,5,8,3,6,9,4,7,10))

# For simulations that included sparse events in the noise structure, PCP-LOD correctly included 69% (IQR: 67%, 71%), 
# 70% (IQR: 68%, 72%), and 65% (IQR: 62%, 67%) of sparse values in the S matrix, on average for simulations with 25%, 
# 50%, and 75% < LOD, respectively.
sim_sparse_id = sim_pcp_out %>% 
  filter(name == "sim_sparse") %>% 
  left_join(sim_lod) %>% 
  mutate(S = map(pcp_out, function(x) x$S)) %>% 
  dplyr::select(-c(true_patterns,true_scores,lod,lod_sqrt2_mat, pcp_out)) %>% 
  mutate(true_sparse_mask = map(sparsity, function(x) (x != 0)),
         lod_mask = map(lod_neg1_mat, function(x) x != -1),
         resid = pmap(list(lod_mask, sim, L), function(lod_mask, sim, L) (sim*lod_mask - L*lod_mask)),
         # do take lod into account bc some sparse events were corrupted to <LOD
         sparse_thresh = map(resid, function(x) apply(x, 2, sd)*2),
         S_thresh_p = map2(S, sparse_thresh, function(S, sparse_thresh) (S > sparse_thresh)*1),
         S_thresh_n = map2(S, sparse_thresh, function(S, sparse_thresh) (S < -sparse_thresh)*-1),
         S_thresh = map2(S_thresh_p, S_thresh_n, function(x,y) abs(x) + abs(y)),
         sparse_prop = map2(true_sparse_mask, S_thresh, function(x,y) sum((x==1) & (y==1))/sum(x==1))) %>% 
  unnest(sparse_prop)

identified_overall = sim_sparse_id  %>% 
  group_by(lim) %>% 
  summarize(quant25 = quantile(sparse_prop, probs = .25), 
            quant50 = quantile(sparse_prop, probs = .5),
            quant75 = quantile(sparse_prop, probs = .75)) %>% 
  mutate_all(round, 2) %>% 
  mutate(chemicals = "Overall") %>% select(chemicals, everything())

identified_size = sim_sparse_id  %>% 
  group_by(chemicals, lim) %>% 
  summarize(quant25 = quantile(sparse_prop, probs = .25), 
            quant50 = quantile(sparse_prop, probs = .5),
            quant75 = quantile(sparse_prop, probs = .75)) %>% 
  mutate_all(round, 2) %>% 
  mutate(chemicals = as.character(chemicals))

sparse_identified = bind_rows(identified_size, identified_overall) %>% 
  pivot_wider(names_from = lim,
              values_from = c(quant25, quant50, quant75)) %>% 
  ungroup() %>% select(c(1, 2,5,8,3,6,9,4,7,10))

# make tables
# sparse_identified %>% 
#   flextable() %>% 
#   add_header_row(colwidths = c(1,3,3,3),
#                  values = c("", "25% < LOD",	"50% < LOD",	"75% < LOD")) %>% 
#   flextable::save_as_docx(path = 'Figures/identified_sparse_table.docx')

# sparse_percent %>% 
#   flextable() %>% 
#   add_header_row(colwidths = c(1,3,3,3),
#                  values = c("", "25% < LOD",	"50% < LOD",	"75% < LOD")) %>% 
#   flextable::save_as_docx(path = 'Figures/percent_sparse_table.docx')
