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
load("./sims/pcp_out_re.rda")

table(unlist(pcp_out_re$pca_rank))
summary(unlist(pcp_out_re$pca_rank))
sum(unlist(pcp_out_re$pca_rank) == 4)/1800

# Get metrics ####
pcp_long = pcp_out_re %>% 
  rename(pcp_pred = L) %>% 
  mutate(pcp_S = map(pcp_out, function(x) x$S),
         pcp_rank = list(4),
         pcp_LS = map2(pcp_pred, pcp_S, `+`),
         pcp_error = map2(chem, pcp_pred, function(x,y) norm(x-y,"F")/norm(x,"F"))) %>% 
  rename_all(~str_replace(., "svd_pcp", "pcp_svd")) %>% 
  rename_all(~str_replace(., "svd_pca", "pca_svd")) %>% 
  rename_all(~str_replace(., "_left", "left")) %>% 
  rename_all(~str_replace(., "_right", "right")) %>% 
  pivot_longer(c(grep("(pca|pcp)", names(.))),
               names_to = c("method", "which"),
               names_sep = "_") %>% 
  pivot_wider(names_from = which,
              values_from = value)

pcp_metrics = pcp_long %>% 
  mutate(mask = map(lod_neg1_mat, function(x) x != -1),
         relerr_above = pmap(list(mask, chem, pred), 
                               function(mask, x,y) norm(mask*x-mask*y,"F")/norm(mask*x,"F")),
         relerr_below= pmap(list(mask, chem, pred), 
                               function(mask, x,y) norm((!mask)*x-(!mask)*y,"F")/norm((!mask)*x,"F")),
         error_right = map2(svd_chemright, svdright, 
                             function(x,y) norm(x-y,"F")/norm(x,"F")),
         error_left = map2(svd_chemleft, svdleft,
                            function(x,y) norm(x-y,"F")/norm(x,"F"))) %>% 
  unnest(c(error, relerr_above, relerr_below, error_right, error_left)) 

summary(unlist(pcp_metrics$relerr_above))
summary(unlist(pcp_metrics$relerr_below))

pcp_metrics %>% 
  # filter(method == "pcp") %>% 
  # pivot_longer(relerr_above:relerr_below,
  #              names_to = "ab") %>% 
  group_by(lim, method, name) %>% 
  summarize(qs = quantile(error), prop = seq(0, 1, 0.25)) %>% 
  pivot_wider(names_from = prop,
              values_from = qs) %>%
  arrange(name, lim) %>% 
  print(n=36)

pcp_metrics %>% 
  pivot_longer(error_left:error_right,
               names_to = "which") %>%
  filter(grepl("left", which)) %>% 
  group_by(chemicals, lim, method, name) %>% 
  summarize(qs = quantile(value), prop = seq(0, 1, 0.25)) %>% 
  pivot_wider(names_from = prop,
              values_from = qs) %>%
  arrange(chemicals, method, name) %>% 
  print(n=36)

# PLOT ####
# Overall error figure
pdf("./sims/sim_boxplots_16.pdf")
pcp_metrics %>%
  filter(chemicals == 16) %>% 
  mutate(lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         method = ifelse(method == "pca", "PCA", "PCP-LOD"),
         name = case_when(name == "sim_1" ~ "N(0, 1)",
                          name == "sim_5" ~ "N(0, 5)",
                          name == "sim_sparse" ~ "N(0, 1) + sparse")) %>% 
  ggplot(aes(x = lim, y = error, color = method, fill = method)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.25, alpha = 0.4) +
  scale_y_log10() +
  facet_grid(.~name) +
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") + 
  scale_color_manual(values = c("#E58601",
                                "#2196F3")) +
  scale_fill_manual(values = c("#E58601",
                               "#2196F3")) +
  theme(legend.title = element_blank())
dev.off()

pdf("./sims/sim_boxplots_48.pdf")
pcp_metrics %>%
  filter(chemicals != 16) %>% 
  mutate(lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         method = ifelse(method == "pca", "PCA", "PCP-LOD"),
         name = case_when(name == "sim_1" ~ "N(0, 1)",
                          name == "sim_5" ~ "N(0, 5)",
                          name == "sim_sparse" ~ "N(0, 1) + sparse")) %>% 
  ggplot(aes(x = lim, y = error, color = method, fill = method)) +
  geom_boxplot(outlier.size = 0.25, alpha = 0.4) +
  scale_y_log10() +
  facet_grid(.~name) +
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") + 
  scale_color_manual(values = c("#E58601",
                                "#2196F3")) +
  scale_fill_manual(values = c("#E58601",
                               "#2196F3")) +
  theme(legend.title = element_blank())
dev.off()

# error above and below lod
pdf("./sims/lod_boxplots_16.pdf", height = 10)
pcp_metrics %>% 
  filter(chemicals == 16) %>% 
  pivot_longer(relerr_above:relerr_below,
               names_to = "which") %>%
  mutate(which = fct_inorder(ifelse(which == "relerr_above", "Above LOD", "Below LOD")),
         lim = case_when(lim == 0.25 ~ "25%",
                                lim == 0.5 ~ "50%",
                                lim == 0.75 ~ "75%"),
                method = ifelse(method == "pca", "PCA", "PCP-LOD"),
                name = case_when(name == "sim_1" ~ "N(0, 1)",
                                 name == "sim_5" ~ "N(0, 5)",
                                 name == "sim_sparse" ~ "N(0, 1) + sparse")) %>% 
  ggplot(aes(x = lim, y = value, color = method, fill = method)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.25, alpha = 0.4) +
  scale_y_log10() +
  facet_grid(which~name) + 
  scale_color_manual(values = c("#E58601", # PCA
                                "#2196F3")) +
  scale_fill_manual(values = c("#E58601", # PCA
                               "#2196F3")) +
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") +
  theme(legend.title = element_blank())
dev.off()

pdf("./sims/lod_boxplots_48.pdf", height = 10)
pcp_metrics %>% 
  filter(chemicals == 48) %>% 
  pivot_longer(relerr_above:relerr_below,
               names_to = "which") %>%
  mutate(which = fct_inorder(ifelse(which == "relerr_above", "Above LOD", "Below LOD")),
         lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         method = ifelse(method == "pca", "PCA", "PCP-LOD"),
         name = case_when(name == "sim_1" ~ "N(0, 1)",
                          name == "sim_5" ~ "N(0, 5)",
                          name == "sim_sparse" ~ "N(0, 1) + sparse")) %>% 
  ggplot(aes(x = lim, y = value, color = method, fill = method)) +
  geom_boxplot(outlier.size = 0.25, alpha = 0.4) +
  scale_y_log10() +
  facet_grid(which~name) + 
  scale_color_manual(values = c("#E58601", # PCA
                                "#2196F3")) +
  scale_fill_manual(values = c("#E58601", # PCA
                               "#2196F3")) +
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") +
  theme(legend.title = element_blank())
dev.off()

# table
table_err = pcp_metrics %>% 
  dplyr::select(seed, which = name, method, lim, error, relerr_above, relerr_below) %>% 
  mutate(lim = str_c(lim, "lim")) %>% 
  group_by(which, lim, method) %>% 
  summarize(qerr = quantile(error), 
            qabove = quantile(relerr_above),
            qbelow = quantile(relerr_below),
            props = seq(0,1,.25)) %>%
  filter(props != 0 & props != 1) %>% 
  pivot_longer(qerr:qbelow) %>% 
  arrange(name) %>% 
  pivot_wider(names_from = c(lim, props),
              values_from = value)
table_err

xtable::xtable(table_err)

# SVD right and left
# pdf("./sims/figures/svd_boxplots_left.pdf", height = 10)
pcp_metrics %>% 
  pivot_longer(error_left:error_right,
               names_to = "which") %>%
  mutate(which = fct_inorder(ifelse(which == "error_left", "Left eigenvectors", "Right eigenvectors")),
         lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         method = ifelse(method == "pca", "PCA", "PCP-LOD"),
         name = case_when(name == "sim_1" ~ "N(0, 1)",
                          name == "sim_5" ~ "N(0, 5)",
                          name == "sim_sparse" ~ "N(0, 1) + sparse"),
         chemicals = str_c(chemicals, " chemicals")) %>%
  filter(grepl("Left", which)) %>% 
  ggplot(aes(x = lim, y = value, color = method, fill = method)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.25, alpha = 0.4) +
  scale_y_log10() +
  facet_grid(chemicals~name) + 
  scale_color_manual(values = c("#E58601", # PCA
                                "#2196F3")) +
  scale_fill_manual(values = c("#E58601", # PCA
                               "#2196F3")) +
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") +
  theme(legend.title = element_blank())
# dev.off()

# pdf("./sims/figures/svd_boxplots_right.pdf", height = 10)
pcp_metrics %>% 
  pivot_longer(error_left:error_right,
               names_to = "which") %>%
  mutate(which = fct_inorder(ifelse(which == "error_left", "Left eigenvectors", "Right eigenvectors")),
         lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         method = ifelse(method == "pca", "PCA", "PCP-LOD"),
         name = case_when(name == "sim_1" ~ "N(0, 1)",
                          name == "sim_5" ~ "N(0, 5)",
                          name == "sim_sparse" ~ "N(0, 1) + sparse"),
         chemicals = str_c(chemicals, " chemicals")) %>%
  filter(grepl("Right", which)) %>% 
  ggplot(aes(x = lim, y = value, color = method, fill = method)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.25, alpha = 0.4) +
  scale_y_log10() +
  facet_grid(chemicals~name, scales = "free_y") + 
  scale_color_manual(values = c("#E58601", # PCA
                                "#2196F3")) +
  scale_fill_manual(values = c("#E58601", # PCA
                               "#2196F3")) +
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") +
  theme(legend.title = element_blank())
# dev.off()

# Sparsity ####
all_sparse = pcp_out_re %>% 
  mutate(S = map(pcp_out, function(x) x$S)) %>% 
  dplyr::select(1:10, L, S) %>% 
  full_join(., sim_lod) %>% 
  mutate(mask = map(lod_neg1_mat, function(x) x != -1), # detected
         resid = pmap(list(mask, lod_neg1_mat, L), function(x,a,b) a*x - b*x),
         sparse_threshold = map(resid, function(x) apply(x,2,sd)*2),
         s_positive = map2(S, sparse_threshold, function(s,t) s > t),
         s_negative = map2(S, sparse_threshold, function(s,t) (s < -t)*-1),
         s_events = map2(s_positive, s_negative, `+`),
         added_sparse_mask = map(sparsity, function(x) x != 0),
         sparse_matrix_mask = map(s_events, abs),
         sparse_in_sparse = map2(added_sparse_mask, s_positive, function(x,y) sum(x*y)/sum(x)),
         percent_sparse = map(sparse_matrix_mask, function(x) sum(x)/(ncol(x)*nrow(x))),
         percent_low = map(s_negative, function(x) sum(x)/(ncol(x)*nrow(x))),
         percent_high = map(s_positive, function(x) sum(x)/(ncol(x)*nrow(x))))
    
any(S0 < 0)

S0 = all_sparse$s_events[[1]]
S0mask = all_sparse$sparse_matrix_mask[[1]]
addedSmask = all_sparse$added_sparse_mask[[1]] * 1

heatmaply::heatmaply(S0)
heatmaply::heatmaply(S0mask)
heatmaply::heatmaply(addedSmask)

heatmaply(S0, show_dendrogram = F, ylab = "Participants", Rowv = F, Colv = F,
          showticklabels = c(T, F),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
            low = "#67A9CF",
            mid= "white",
            high = "#EF8A62", 
            midpoint = 0))

heatmaply(addedSmask, show_dendrogram = F, ylab = "Participants", Rowv = F, Colv = F,
          showticklabels = c(T, F),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
            low = "#67A9CF",
            mid= "white",
            high = "#EF8A62", 
            midpoint = 0))

all_sparse %>% 
  filter(name == "sim_sparse") %>% 
  unnest(sparse_in_sparse) %>% 
  group_by(lim) %>% 
  summarize(qs = quantile(sparse_in_sparse), probs = seq(0,1,.25)) %>% 
  pivot_wider(names_from = probs,
              values_from = qs)

all_sparse %>% 
  unnest(percent_sparse) %>% 
  group_by(lim) %>% 
  summarize(qs = mean(percent_sparse)) %>% 
  #summarize(qs = quantile(percent_sparse), probs = seq(0,1,.25)) %>% 
  pivot_wider(names_from = probs,
              values_from = qs)

sum(S0mask)/(500*16)
sum(addedSmask)/(500*16)
sum(addedSmask*S0mask)/(500*16)

sum(addedSmask*S0mask)/sum(addedSmask)

sum(abs(Smask * S0) > 0.0001)/(500*16)



