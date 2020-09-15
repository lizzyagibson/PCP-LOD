# Plot PCP-LOD, LOD/sqrt(2), PCA resulting error
# Across rank, sigma, mu, and LOD
library(tidyverse)
library(R.matlab)
library(RColorBrewer)
library(wesanderson)
theme_set(theme_bw(base_size = 20) + theme(legend.position = "bottom",
                             strip.background =element_rect(fill="white")))
rank = c(1, 2, 3, 4, 5)
sigma = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, "1")
delta = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, "1")

length(rank)*length(sigma)*length(delta)

out_mu <- tibble()
pcp_out <- tibble()

##### Read loop
for (r in rank) {
for (d in delta) {
for (s in sigma) {
  
  if (file.exists(paste0("HPC_PCP/mat_out/pcp_all_lod_", d, "_rank_", r, 
                                    "_sigma_", s, ".mat"))) {
  read <- readMat(here::here(paste0("HPC_PCP/mat_out/pcp_all_lod_", d, "_rank_", r, 
                                    "_sigma_", s, ".mat")))}
  
  relerror <- read[[1]] %>% 
    as_tibble() %>% rename("relerror_PCA" =1, "relerror_PCP_sqrt2" =2, "relerror_PCP_LOD" =3)

  out_mu <- read[[2]] %>% 
    as_tibble() %>% rename("belowlod_relerror_PCA" =1, "belowlod_relerror_PCP_sqrt2" =2, 
                           "belowlod_relerror_PCP_LOD" =3) %>% 
    mutate(rank = as.factor(r),
           sigma = as.factor(s),
           delta = as.factor(d),
           iter = 1:100) %>% 
    cbind(., relerror) %>% as_tibble()
  
  pcp_out <- rbind(pcp_out, out_mu) 
  }
}
}
#####

pcp_plot <- pcp_out %>%
  select(rank, sigma, delta, iter, belowlod_relerror_PCA:belowlod_relerror_PCP_LOD, 
         relerror_PCA:relerror_PCP_LOD) %>% 
  pivot_longer(grep("relerror", colnames(.)),
               values_to = "rel_error") %>% 
  mutate(type = ifelse(grepl("belowlod", name), "belowlod", 
                       ifelse(grepl("score", name), "score",
                       ifelse(grepl("pattern", name), "pattern", "overall"))),
         model = ifelse(grepl("PCA", name), "PCA",
                        ifelse(grepl("sqrt2", name), "PCP_sqrt2", "PCP_LOD")),
         rank = as.factor(rank),
         sigma = as.factor(sigma),
         delta = as.factor(delta)) %>% 
  select(-name) %>% 
  mutate(type = case_when(type == "belowlod" ~ "< LOD",
                          type == "overall" ~ "Overall"),
         type = fct_relevel(type, "Overall", "< LOD")) %>% 
  mutate(model = case_when(model == "PCA" ~ "PCA",
                           model == "PCP_LOD" ~ "PCP-LOD",
                           model == "PCP_sqrt2" ~ #paste0("PCP w/ LOD/", intToUtf8(8730),2))) %>% 
                             "PCP"))
pcp_plot

pcp_plot %>% 
  filter(delta != "1") %>% 
  filter(sigma == "0") %>% 
  filter(rank == "1") %>%
  ggplot(aes(x = delta, y = rel_error, color = model)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.25) +
  facet_grid(~type) + scale_y_log10() +
  scale_color_manual(values = c("#E58601", # PCA
                                "#B40F20",  # PCP SQRT2
                                "#2196F3" # PCP LOD
  )) +
  #scale_color_manual(values = wes_palette(n=3, name="Darjeeling1")) +
  labs(x = "Values < LOD (proportion)",
       y = "Relative Prediction Error") +
  theme(
    legend.title = element_blank())





