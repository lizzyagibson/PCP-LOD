#####
#####
# Plot PCP-LOD, LOD/sqrt(2), PCA resulting error
# Across rank, sigma, mu, and LOD
# 9/16/2020
# PCP-LOD with adaptive rho
#####
#####

## Load libraries
library(tidyverse)
library(R.matlab)
library(RColorBrewer)
library(ggsci)

## Set plot theme
theme_set(theme_bw(base_size = 20) + 
            theme(legend.position = "bottom",
                             strip.background =element_rect(fill="white"),
              legend.title = element_blank()))

## Loop these
rank = c(1, 2, 3, 4, 5)
sigma = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, "1")
delta = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, "1")

## Array length for HPC
length(rank)*length(sigma)*length(delta)

out_mu <- tibble()
pcp_out <- tibble()

## Read loop
for (r in rank) {
for (d in delta) {
for (s in sigma) {
  
  if (file.exists((here::here(paste0("HPC_PCP/mat_out/pcp_all_lod_", d, "_rank_", r, 
                                    "_sigma_", s, ".mat"))))) {
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

pcp_out

## Arrange data for ggplot
pcp_plot <- pcp_out %>%
  pivot_longer(grep("relerror", colnames(.)),
               values_to = "rel_error") %>% 
  mutate(type = ifelse(grepl("belowlod", name), "< LOD", "Overall"),
         model = case_when(grepl("PCA", name) ~ "PCA",
                           grepl("sqrt2", name) ~ "PCP_sqrt2", 
                           TRUE ~ "PCP_LOD")) %>% 
  mutate_at(vars(1:3), as.factor) %>% 
  select(-name) %>% 
  mutate(type = fct_relevel(type, "Overall", "< LOD"))

pcp_plot

pcp_plot %>% 
  filter(!(delta %in% c("0.7", "0.8", "0.9", "1"))) %>% # dont really care
  filter(sigma == "0") %>% # Change
  filter(rank == "1") %>% # Change
  ggplot(aes(x = delta, y = rel_error, color = model)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.25) +
  facet_grid(~type) + scale_y_log10() +
  scale_color_nejm() +
  labs(x = "Values < LOD (proportion)",
       y = "Relative Prediction Error")





