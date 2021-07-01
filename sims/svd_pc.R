# Get SVD of PCP solution

# Load functions
source("Sims/functions.R")

# Load PCP & PCA results
load("./sims/pcp_out_all.rda")

# SVD ####
svd_out = pcp_out %>% 
  mutate(svd_chem = map(chem, svd),
         svd_pca  = map(lod_sqrt2_mat, svd),
         L        = map(pcp_out, function(x) x$L),
         svd_pcp  = map(L, svd),
         svd_chem_left  = map(svd_chem, function(x) x$u),
         svd_chem_right = map(svd_chem, function(x) x$v),
         svd_pca_left  = map(svd_pca, function(x) x$u),
         svd_pca_right = map(svd_pca, function(x) x$v),
         svd_pcp_left  = map(svd_pcp, function(x) x$u),
         svd_pcp_right = map(svd_pcp, function(x) x$v))

# Rearrange vectors to be closest to original SVD results
svd_re = svd_out %>% 
  mutate(svd_pca_left  = map2(svd_chem_left,  svd_pca_left,  factor_correspondence),
         svd_pcp_left  = map2(svd_chem_left,  svd_pcp_left,  factor_correspondence),
         svd_pcp_right = map2(svd_chem_right, svd_pcp_right, factor_correspondence),
         svd_pca_right = map2(svd_chem_right, svd_pca_right, factor_correspondence))

# save(svd_re, file = "./Sims/Sim Data/svd_out_re.rda")  
load("./Sims/Sim Data/svd_out_re.rda")  

svd_metrics = svd_re %>% 
  mutate(svd_err_pca_left   = map2(svd_chem_left, svd_pca_left, function(x,y) norm(x-y,"F")/norm(x,"F")),
         svd_err_pcp_left   = map2(svd_chem_left, svd_pcp_left, function(x,y) norm(x-y,"F")/norm(x,"F")),
         svd_err_pca_right   = map2(svd_chem_right, svd_pca_right, function(x,y) norm(x-y,"F")/norm(x,"F")),
         svd_err_pcp_right   = map2(svd_chem_right, svd_pcp_right, function(x,y) norm(x-y,"F")/norm(x,"F"))) %>% 
  unnest(c(svd_err_pca_left, svd_err_pcp_left, svd_err_pca_right, svd_err_pcp_right)) %>% 
  group_by(chemicals, lim) %>% 
  summarize(pca_left = quantile(svd_err_pca_left), prop = seq(0,1,.25),
            pcp_left = quantile(svd_err_pcp_left),
            pca_right = quantile(svd_err_pca_right),
            pcp_right = quantile(svd_err_pcp_right))

svd_metrics %>% 
  pivot_longer(c(pca_left, pcp_left:pcp_right)) %>% 
  pivot_wider(names_from = prop,
              values_from = value) %>% print(n=24)
