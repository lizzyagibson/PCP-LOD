
# Load experiment results
load("./Sims/Sim Data/pca_metrics.rda")
load("./Sims/Sim Data/pca_svd_metrics.rda")

load("./Sims/Sim Data/pcp_metrics.rda")
load("./Sims/Sim Data/pcp_svd_metrics.rda")

# Load functions
source("./functions.R")

# Combine
metrics     = full_join(pcp_metrics, pca_metrics)
svd_metrics = full_join(pcp_svd_metrics, pca_svd_metrics)

# PLOT ####

# make colorblind friendly!
display.brewer.all(colorblindFriendly = TRUE)

# Overall error figure

# Figure 2 in manuscript
#pdf("./Figures/Figure_2.pdf")
metrics %>%
  filter(chemicals == 16) %>% 
  mutate(lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         name = case_when(name == "sim_1" ~ "Low Noise",
                          name == "sim_5" ~ "High Noise",
                          name == "sim_sparse" ~ "Sparse Events"),
         name = fct_relevel(name, "Low Noise", "Sparse Events", "High Noise")) %>% 
  # Jaime: in the manuscript (figure 2?) facets are named differently, please make them coincide
  # Lizzy: Thanks, updated manuscript to have these labels
  ggplot(aes(x = lim, y = all_relerr, color = method, fill = method)) +
  geom_boxplot_pattern(aes(pattern=method, pattern_fill = method, 
                           pattern_colour = method, pattern_density = method),
                       notch = TRUE, outlier.size = 0.25, alpha = 0.4,
                       pattern_spacing = 0.015) +
  scale_y_log10() +
  facet_grid(.~name) +
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") + 
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  scale_pattern_fill_brewer(palette = "Dark2") +
  scale_pattern_color_brewer(palette = "Dark2") +
  scale_pattern_density_manual(values = c(0.1, 0)) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5, 'cm'))
#dev.off()

#pdf("./Figures/Figure_S1.pdf")
# Jaime: named as Figure S1 in manuscript
metrics %>%
  filter(chemicals != 16) %>% 
  mutate(lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         name = case_when(name == "sim_1" ~ "Low Noise",
                          name == "sim_5" ~ "High Noise",
                          name == "sim_sparse" ~ "Sparse Events"),
         name = fct_relevel(name, "Low Noise", "Sparse Events", "High Noise")) %>% 
  ggplot(aes(x = lim, y = all_relerr, color = method, fill = method)) +
  geom_boxplot_pattern(aes(pattern=method, pattern_fill = method,
                           pattern_colour = method, pattern_density = method),
                       outlier.size = 0.25, alpha = 0.4,
                       pattern_spacing = 0.015) +
  scale_y_log10() +
  facet_grid(.~name) +
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") + 
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  scale_pattern_fill_brewer(palette = "Dark2") +
  scale_pattern_color_brewer(palette = "Dark2") +
  scale_pattern_density_manual(values = c(0.1, 0)) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5, 'cm'))
#dev.off()

# Figure 3 in manuscript
# error above and below lod
#pdf("./Figures/Figure_3.pdf", height = 10)
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
  # Jaime: in the manuscript (figure 3?) facets are named differently, please make them coincide
  # Lizzy: Thanks, I updated the manuscript to match these
  ggplot(aes(x = lim, y = value, color = method, fill = method)) +
  geom_boxplot_pattern(aes(pattern=method, pattern_fill = method,
                           pattern_colour = method, pattern_density = method),
                       notch = TRUE, outlier.size = 0.25, alpha = 0.4,
                       pattern_spacing = 0.015) +
  scale_y_log10() +
  facet_grid(which~name) + 
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  scale_pattern_fill_brewer(palette = "Dark2") +
  scale_pattern_color_brewer(palette = "Dark2") +
  scale_pattern_density_manual(values = c(0.1, 0)) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5, 'cm'))
#dev.off()

#pdf("./Figures/Figure_S2.pdf", height = 10)
# Jaime: named as Figure S2 in manuscript
metrics %>% 
  filter(chemicals == 48) %>% 
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
  geom_boxplot_pattern(aes(pattern=method, pattern_fill = method,
                           pattern_colour = method, pattern_density = method),
                       outlier.size = 0.25, alpha = 0.4,
                       pattern_spacing = 0.015) +
  scale_y_log10() +
  facet_grid(which~name) + 
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  scale_pattern_fill_brewer(palette = "Dark2") +
  scale_pattern_color_brewer(palette = "Dark2") +
  scale_pattern_density_manual(values = c(0.1, 0)) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5, 'cm'))
#dev.off()

# tables ####
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
# Jaime: named as Table S1 in manuscript
xtable::xtable(table_err)

# Corresponding numeric data for figures 2 & S1
# Relative predictive error of PCP-LOD and PCA stratified by p (16 or 48)
table_fig2_s1 = metrics %>% 
  dplyr::select(seed, which_sim = name, method, lim, chemicals, all_relerr) %>% 
  mutate(lim = str_c(lim, "lim")) %>% 
  group_by(which_sim, lim, method, chemicals) %>% 
  summarize(qerr = quantile(all_relerr, na.rm = T), 
            props = c("Min", "Q25", "Median", "Q75", "Max")) %>%
  #filter(props != 0 & props != 1) %>% 
  pivot_longer(qerr) %>% 
  arrange(name) %>% 
  pivot_wider(names_from = c(lim, props),
              values_from = value) %>% 
  select(-name) %>% 
  select(which_sim, chemicals, method, everything()) %>% 
  arrange(which_sim, chemicals, method) %>% 
  mutate_if(is.numeric, round, 2)

# table_fig2_s1 %>% flextable() %>% flextable::save_as_docx(path = 'Tables/table_fig2_s1.docx')

# Corresponding numeric data for figures 3 & S1
# Relative predictive error of PCP-LOD and PCA stratified by p (16 or 48) and detection
table_fig3_s2 = metrics %>% 
  dplyr::select(seed, which_sim = name, method, lim, chemicals, above, below) %>% 
  mutate(lim = str_c(lim, "lim")) %>% 
  group_by(which_sim, lim, method, chemicals) %>% 
  summarize(qabove = quantile(above, na.rm = T), 
            qbelow = quantile(below, na.rm = T), 
            props = c("Min", "Q25", "Median", "Q75", "Max")) %>%
  #filter(props != 0 & props != 1) %>% 
  pivot_longer(qabove:qbelow) %>% 
  arrange(name) %>% 
  pivot_wider(names_from = c(lim, props),
              values_from = value) %>% 
  select(name, which_sim, chemicals, method, everything()) %>% 
  arrange(name, which_sim, chemicals, method) %>% 
  mutate_if(is.numeric, round, 2)

# table_fig3_s2 %>% flextable() %>% flextable::save_as_docx(path = 'Tables/table_fig3_s2.docx')

# Corresponding numeric data for figures 3 & S2
# Relative predictive error of PCP-LOD and PCA stratified by detection and p (16 or 48)
table_fig4_s3 = svd_metrics %>% 
  dplyr::select(seed, which_sim = name, method, lim, chemicals, left, right) %>% 
  mutate(lim = str_c(lim, "lim")) %>% 
  group_by(which_sim, lim, method, chemicals) %>% 
  summarize(qleft = quantile(left, na.rm = T), 
            qright = quantile(right, na.rm = T), 
            props = c("Min", "Q25", "Median", "Q75", "Max")) %>%
  #filter(props != 0 & props != 1) %>% 
  pivot_longer(qleft:qright) %>% 
  arrange(name) %>% 
  pivot_wider(names_from = c(lim, props),
              values_from = value) %>% 
  select(name, which_sim, chemicals, method, everything()) %>% 
  arrange(name, which_sim, chemicals, method) %>% 
  mutate_if(is.numeric, round, 2)

# table_fig4_s3 %>% flextable() %>% flextable::save_as_docx(path = 'Tables/table_fig4_s3.docx')

# more figures ####

# SVD right and left
#pdf("./Figures/Figure_4.pdf", height = 10)
# Jaime: named as Figure 4 in manuscript
svd_metrics %>% 
  pivot_longer(left:right,
               names_to = "which") %>%
  mutate(which = fct_inorder(ifelse(which == "left", "Left eigenvectors", "Right eigenvectors")),
         lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         chemicals = str_c(chemicals, " chemicals"),
         name = case_when(name == "sim_1" ~ "Low Noise",
                          name == "sim_5" ~ "High Noise",
                          name == "sim_sparse" ~ "Sparse Events"),
         name = fct_relevel(name, "Low Noise", "Sparse Events", "High Noise")) %>% 
  filter(grepl("Left", which)) %>% 
  ggplot(aes(x = lim, y = value, color = method, fill = method)) +
  geom_boxplot_pattern(aes(pattern=method, pattern_fill = method,
                           pattern_colour = method, pattern_density = method),
                       notch = TRUE, outlier.size = 0.25, alpha = 0.4,
                       pattern_spacing = 0.015) +
  scale_y_log10() +
  facet_grid(chemicals~name) + 
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  scale_pattern_fill_brewer(palette = "Dark2") +
  scale_pattern_color_brewer(palette = "Dark2") +
  scale_pattern_density_manual(values = c(0.1, 0)) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5, 'cm'))
#dev.off()

#pdf("./Figures/Figure_S3.pdf", height = 10)
# Jaime: named as Figure S3 in manuscript
svd_metrics %>% 
  pivot_longer(left:right,
               names_to = "which") %>%
  mutate(which = fct_inorder(ifelse(which == "error_left", "Left eigenvectors", "Right eigenvectors")),
         lim = case_when(lim == 0.25 ~ "25%",
                         lim == 0.5 ~ "50%",
                         lim == 0.75 ~ "75%"),
         chemicals = str_c(chemicals, " chemicals"),
         name = case_when(name == "sim_1" ~ "Low Noise",
                          name == "sim_5" ~ "High Noise",
                          name == "sim_sparse" ~ "Sparse Events"),
         name = fct_relevel(name, "Low Noise", "Sparse Events", "High Noise")) %>% 
  filter(grepl("Right", which)) %>% 
  ggplot(aes(x = lim, y = value, color = method, fill = method)) +
  geom_boxplot_pattern(aes(pattern=method, pattern_fill = method,
                           pattern_colour = method, pattern_density = method),
                       outlier.size = 0.25, alpha = 0.4,
                       pattern_spacing = 0.015) +
  scale_y_log10() +
  facet_grid(chemicals~name, scales = "free_y") + 
  labs(x = "Percentage of values < LOD",
       y = "Relative Prediction Error") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  scale_pattern_fill_brewer(palette = "Dark2") +
  scale_pattern_color_brewer(palette = "Dark2") +
  scale_pattern_density_manual(values = c(0.1, 0)) +
  theme(legend.title = element_blank(),
        legend.key.size = unit(1.5, 'cm'))
#dev.off()
