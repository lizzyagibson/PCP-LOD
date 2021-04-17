# CV

library(tidyverse)
library(GGally)
library(PCPhelpers)
library(pcpr)
library(patchwork)
library(factoextra)
library(CVXR)

theme_set(theme_bw(base_size = 20) + theme(legend.position = "bottom",
                                           strip.background =element_rect(fill="white")))

# 4 patterns
# 2 mixture sizes, 16 & 48
# 1 sample size = 500
# 3 proportions <LOD, 25, 50, 75

load("./sims/sim_lod.RDA")
sim_lod

# Run PCP-LOD ####
# Run PCP-LOD on all sims!
# sim_pcp_out = sim_lod %>% 
#               mutate(pcp_out = map2(lod_neg1_mat, lod, function(x,y)
#                         root_pcp_noncvx_nonnegL_na_lod(D = x, MAX_ITER = 20000,
#                                                   lambda = 1/sqrt(nrow(x)), mu = sqrt(ncol(x)/2), 
#                                                   r = 4, LOD = y,
#                                                           verbose = TRUE)))
# 
# pcp_out = tibble()
# for (i in 1:600) {
#   load(paste0("/ifs/scratch/msph/ehs/eag2186/pcp/lod/sim_pcp_out_", i, ".rda"))
#   pcp_out <- bind_rows(pcp_out, sim_pcp_out)
#        print(i)
# }
load("./sims/pcp_out_all.rda")

# SVD ####
svd_out = pcp_out %>% 
  mutate(svd_chem = map(chem, svd),
         svd_pca  = map(lod_sqrt2_mat, svd),
         L = map(pcp_out, function(x) x$L),
         svd_pcp  = map(L, svd),
         svd_chem_left  = map(svd_chem, function(x) x$u),
         svd_chem_right = map(svd_chem, function(x) x$v),
         svd_pca_left  = map(svd_pca, function(x) x$u),
         svd_pca_right = map(svd_pca, function(x) x$v),
         svd_pcp_left  = map(svd_pcp, function(x) x$u),
         svd_pcp_right = map(svd_pcp, function(x) x$v),
         svd_pca_left  = map2(svd_chem_left, svd_pca_left, factor_correspondence),
         svd_pca_right = map2(svd_chem_right, svd_pca_right, factor_correspondence),
         svd_pcp_left  = map2(svd_chem_left, svd_pcp_left, factor_correspondence),
         svd_pcp_right = map2(svd_chem_right, svd_pcp_right, factor_correspondence))
save(svd_out, file = "./sims/svd_out.rda")  

# Get metrics ####
pcp_metrics = pcp_out %>% 
  mutate(method = "PCP-LOD") %>% 
        arrange(seed, chem) %>% 
            mutate(L = map(pcp_out, function(x) x$L),
                   S = map(pcp_out, function(x) x$S),
                   LS = map2(L,S,`+`)) %>% 
  mutate(relerr_l_all = map2(chem, L, function(x,y) norm(x-y,"F")/norm(x,"F")),
         relerr_ls_all   = map2(chem, LS, function(x,y) norm(x-y,"F")/norm(x,"F"))) %>% 
  unnest(c(relerr_l_all, relerr_ls_all)) %>% 
  mutate(mask = map(lod_neg1_mat, function(x) x != -1),
         relerr_l_above = pmap(list(mask, chem, L), 
                               function(mask, x,y) norm(mask*x-mask*y,"F")/norm(mask*x,"F")),
         relerr_ls_above = pmap(list(mask, chem, LS), 
                               function(mask, x,y) norm(mask*x-mask*y,"F")/norm(mask*x,"F")),
         relerr_l_below = pmap(list(mask, chem, L), 
                               function(mask, x,y) norm((!mask)*x-(!mask)*y,"F")/norm((!mask)*x,"F")),
         relerr_ls_below= pmap(list(mask, chem, L), 
                               function(mask, x,y) norm((!mask)*x-(!mask)*y,"F")/norm((!mask)*x,"F"))) %>% 
  unnest(c(relerr_l_above, relerr_ls_above, relerr_l_below, relerr_ls_below)) 

summary(pcp_metrics$relerr_ls_all)
summary(pcp_metrics$relerr_l_all)

# Combine ####
combined = pcp_metrics %>% 
  dplyr::select(seed, chemicals, lim, method, relerr_l_all, relerr_ls_all, 
                relerr_l_above, relerr_ls_above, relerr_l_below, relerr_ls_below) %>% 
  pivot_longer(relerr_l_all:relerr_ls_below,
               names_to = c("metric", "LS", "which"),
               names_sep = "_") %>% 
  pivot_wider(names_from = c(which,metric),
              values_from = value) %>% 
  mutate(LS = ifelse(LS == "l", "PCP-LOD L", "PCP-LOD L+S")) %>% 
  bind_rows(., pca_combine)

combined %>% 
  group_by(chemicals, lim, method, LS) %>% 
  summarize(qs = quantile(all_relerr), prop = seq(0, 1, 0.25)) %>% 
  pivot_wider(names_from = prop,
              values_from = qs)

# Overall error figure
combined %>% 
  mutate_at(vars(2:5), as_factor) %>%
  mutate(lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%")) %>% 
  ggplot(aes(x = lim, y = all_relerr, color = LS, fill = LS)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.25, alpha = 0.4) +
  scale_y_log10() +
  #facet_wrap(.~chemicals) +
  scale_color_manual(values = c("#E58601", # PCA
                                "#B40F20",  # PCP SQRT2
                                "#2196F3" # PCP LOD
  )) +
  scale_fill_manual(values = c("#E58601", # PCA
                                "#B40F20",  # PCP SQRT2
                                "#2196F3" # PCP LOD
  )) +
  #scale_color_manual(values = wes_palette(n=3, name="Darjeeling1")) +
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") +
  theme(
    legend.title = element_blank())

# error above and below lod
combined %>% 
  mutate_at(vars(2:5), as_factor) %>% 
  pivot_longer(above_relerr:below_relerr,
               names_to = c("which", "drop"),
               names_sep = "_") %>% 
  mutate(which = fct_inorder(ifelse(which == "above", "> LOD", "< LOD")),
         lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%")) %>% 
  ggplot(aes(x = lim, y = value, color = LS, fill = LS)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.25, alpha = 0.4) +
  scale_y_log10() +
  facet_grid(.~which) + 
  scale_color_manual(values = c("#E58601", # PCA
                                "#B40F20",  # PCP SQRT2
                                "#2196F3" # PCP LOD
  )) +
  scale_fill_manual(values = c("#E58601", # PCA
                               "#B40F20",  # PCP SQRT2
                               "#2196F3" # PCP LOD
  )) +
  #scale_color_manual(values = wes_palette(n=3, name="Darjeeling1")) +
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") +
  theme(
    legend.title = element_blank())


