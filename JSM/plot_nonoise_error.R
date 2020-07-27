# Plot PCP-LOD, LOD/sqrt(2), PCA resulting error
# Across rank, sigma, mu, and LOD
library(tidyverse)
library(R.matlab)
library(RColorBrewer)
library(wesanderson)
theme_set(theme_minimal() + theme(legend.position = "bottom"))

rank = c(1, 2, 3, 4, 5)
delta = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

# r = 1
# d = 0.3

nonoise_out <- tibble()

for (r in rank) {
  for (d in delta) {
      
      read <- readMat(here::here(paste0("HPC_PCP/nonoise/no_noise_lod_", d, "_rank_", r, ".mat")))
      
      relerror <- read[[1]] %>% 
        as_tibble() %>% rename("relerror_PCP_LOD" =1, "relerror_PCP_sqrt2" =2)
      
      relerror_safe <- read[[3]] %>% 
        as_tibble() %>% rename("safe_relerror_PCP_LOD" =1, "safe_relerror_PCP_sqrt2" =2)
      
      out <- read[[2]] %>% 
        as_tibble() %>% rename("belowlod_relerror_PCP_LOD" =1, "belowlod_relerror_PCP_sqrt2" =2) %>% 
        mutate(rank = r,
               delta = d,
               iter = 1:100) %>% 
        cbind(., relerror, relerror_safe) %>% as_tibble()
      
      
      nonoise_out <- rbind(nonoise_out, out) %>% drop_na(.)
    }
  }


nonoise_plot <- 
  nonoise_out %>% 
  pivot_longer(grep("relerror", colnames(.)),
               values_to = "rel_error") %>% 
  mutate(type = ifelse(grepl("belowlod", name), "belowlod", 
                       ifelse(grepl("safe", name), "safe", "overall")),
         model = ifelse(grepl("PCA", name), "PCA",
                        ifelse(grepl("sqrt2", name), "PCP_sqrt2", "PCP_LOD")),
         rank = as.factor(rank),
         delta = as.factor(delta)) %>% 
  select(-name)

names(nonoise_plot)

# rank increase -- PCP_LOD does better
nonoise_plot %>% 
  #filter(rank == "5") %>% 
  #filter(type != "overall") %>% 
  #filter(type != "belowlod") %>% 
  ggplot(aes(x = delta, y = rel_error, color = model)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  facet_grid(rank~type, scales = "free_y") + scale_y_log10() +
  theme_bw() + scale_color_manual(values=wes_palette(n=3, name="FantasticFox1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

prop <- nonoise_out %>% 
  mutate(Overall = ifelse(relerror_PCP_LOD <= relerror_PCP_sqrt2, 1, 0),
         BelowLOD = ifelse(belowlod_relerror_PCP_LOD <= belowlod_relerror_PCP_sqrt2, 1, 0),
         SafeRows = ifelse(safe_relerror_PCP_LOD <= safe_relerror_PCP_sqrt2, 1, 0)) %>%
  select(rank, delta, iter, Overall, BelowLOD, SafeRows) %>% 
  pivot_longer(Overall:SafeRows) %>% 
  group_by(rank, delta, name) %>%
  summarize(Better = sum(value)/n()) %>% ungroup() %>% 
  mutate(rank = as.factor(rank),
         delta = as.factor(delta))

prop %>% ggplot(aes(y = Better, x = delta)) +
  geom_line(aes(group = rank, color = rank)) +
  facet_wrap(~name) + scale_color_manual(values=wes_palette(n=5, name="FantasticFox1"))




