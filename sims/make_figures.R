# Load experiment results
load("./Sims/Sim Data/pca_metrics.rda")
load("./Sims/Sim Data/pca_svd_metrics.rda")

load("./Sims/Sim Data/pcp_metrics.rda")
load("./Sims/Sim Data/pcp_svd_metrics.rda")

# Load functions
source("./functions.R")

# Combine
metrics = full_join(pcp_metrics, pca_metrics)
svd_metrics = full_join(pcp_svd_metrics, pca_svd_metrics)

# PLOT ####

# Overall error figure

#pdf("./sims/sim_boxplots_16.pdf")
metrics %>%
  filter(chemicals == 16) %>% 
  mutate(lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         name = case_when(name == "sim_1" ~ "Low Noise",
                          name == "sim_5" ~ "High Noise",
                          name == "sim_sparse" ~ "Sparse Events"),
         name = fct_relevel(name, "Low Noise", "Sparse Events", "High Noise")) %>% 
  ggplot(aes(x = lim, y = all_relerr, color = method, fill = method)) +
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
#dev.off()

#pdf("./sims/sim_boxplots_48.pdf")
metrics %>%
  filter(chemicals != 16) %>% 
  mutate(lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         name = case_when(name == "sim_1" ~ "N(0, 1)",
                          name == "sim_5" ~ "N(0, 5)",
                          name == "sim_sparse" ~ "N(0, 1) + sparse")) %>% 
  ggplot(aes(x = lim, y = all_relerr, color = method, fill = method)) +
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
#dev.off()

# error above and below lod
#pdf("./sims/lod_boxplots_16.pdf", height = 10)
metrics %>% 
  filter(chemicals == 16) %>% 
  pivot_longer(above:below,
               names_to = "which") %>%
  mutate(which = fct_inorder(ifelse(which == "above", "Above LOD", "Below LOD")),
         lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         name = case_when(name == "sim_1" ~ "Low Noise",
                          name == "sim_5" ~ "High Noise",
                          name == "sim_sparse" ~ "Sparse Events"),
         name = fct_relevel(name, "Low Noise", "Sparse Events", "High Noise")) %>% 
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
#dev.off()

#pdf("./sims/lod_boxplots_48.pdf", height = 10)
metrics %>% 
  filter(chemicals == 48) %>% 
  pivot_longer(above:below,
               names_to = "which") %>%
  mutate(which = fct_inorder(ifelse(which == "above", "Above LOD", "Below LOD")),
         lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
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
#dev.off()

# table
table_err = metrics %>% 
  dplyr::select(seed, which = name, method, lim, all_relerr, above, below) %>% 
  mutate(lim = str_c(lim, "lim")) %>% 
  group_by(which, lim, method) %>% 
  summarize(qerr = quantile(all_relerr, na.rm = T), 
            qabove = quantile(above, na.rm = T),
            qbelow = quantile(below, na.rm = T),
            props = seq(0,1,.25)) %>%
  filter(props != 0 & props != 1) %>% 
  pivot_longer(qerr:qbelow) %>% 
  arrange(name) %>% 
  pivot_wider(names_from = c(lim, props),
              values_from = value)

table_err
# Latex version
xtable::xtable(table_err)

# SVD right and left
# pdf("./sims/figures/svd_boxplots_left.pdf", height = 10)
svd_metrics %>% 
  pivot_longer(left:right,
               names_to = "which") %>%
  mutate(which = fct_inorder(ifelse(which == "left", "Left eigenvectors", "Right eigenvectors")),
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
svd_metrics %>% 
  pivot_longer(left:right,
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

#pdf("./sims/figures/svd_boxplots_16.pdf", height = 10)
svd_metrics %>% 
  pivot_longer(error_left:error_right,
               names_to = "which") %>%
  mutate(which = fct_inorder(ifelse(which == "error_left", "Individual scores", "Chemical loadings")),
         lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         method = ifelse(method == "pca", "PCA", "PCP-LOD"),
         name = case_when(name == "sim_1" ~ "N(0, 1)",
                          name == "sim_5" ~ "N(0, 5)",
                          name == "sim_sparse" ~ "N(0, 1) + sparse"),
         chemicals = str_c(chemicals, " chemicals")) %>%
  #filter(grepl("Right", which)) %>% 
  filter(chemicals == "16 chemicals") %>% group_by(lim, method, which, name) %>% 
  summarise(m=median(value)) %>% print(n=36)

ggplot(aes(x = lim, y = value, color = method, fill = method)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.25, alpha = 0.4) +
  scale_y_log10() +
  facet_grid(which~name, scales = "free_y") + 
  scale_color_manual(values = c("#E58601", # PCA
                                "#2196F3")) +
  scale_fill_manual(values = c("#E58601", # PCA
                               "#2196F3")) +
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") +
  theme(legend.title = element_blank())
# dev.off()
