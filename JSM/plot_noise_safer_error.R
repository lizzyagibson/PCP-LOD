# Plot PCP-LOD, LOD/sqrt(2), PCA resulting error
# Across rank, sigma, mu, and LOD
library(tidyverse)
library(R.matlab)
library(RColorBrewer)
library(wesanderson)
theme_set(theme_bw(base_size = 20) + theme(legend.position = "bottom",
                                           strip.background =element_rect(fill="white")))

delta = c("0", 0.05, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

r = 4

noise_out <- tibble()
for (d in delta) {
    read <- readMat(here::here(paste0("HPC_PCP/prop_noise/safe_sigma_0.6_lod_", d, "_rank_4.mat")))
    relerror <- read[[1]] %>% 
      as_tibble() %>% rename("relerror_PCP_LOD" =1, "relerror_PCP_sqrt2" =2)
    relerror_safe <- read[[3]] %>% 
      as_tibble() %>% rename("safe_relerror_PCP_LOD" =1, "safe_relerror_PCP_sqrt2" =2)
    out <- read[[2]] %>% 
      as_tibble() %>% rename("belowlod_relerror_PCP_LOD" =1, "belowlod_relerror_PCP_sqrt2" =2) %>% 
      mutate(rank = 4,
             delta = d,
             iter = 1:100) %>% 
      cbind(., relerror, relerror_safe) %>% as_tibble()
    noise_out <- rbind(noise_out, out) %>% drop_na(iter)
}

noise_plot <- 
  noise_out %>% 
  pivot_longer(grep("relerror", colnames(.)),
               values_to = "rel_error") %>% 
  mutate(type = ifelse(grepl("belowlod", name), "belowlod", 
                       ifelse(grepl("safe", name), "safe", "overall")),
         model = ifelse(grepl("PCA", name), "PCA",
                        ifelse(grepl("sqrt2", name), "PCP_sqrt2", "PCP_LOD")),
         rank = as.factor(rank)) %>% 
  select(-name)

names(noise_plot)

noise_plot %>%  
  #filter(type != "overall") %>% 
  #filter(type != "belowlod") %>% 
  ggplot(aes(x = delta, y = rel_error, color = model)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  facet_grid(rank~type, scales = "free_y") + scale_y_log10() +
  scale_color_manual(values=wes_palette(n=3, name="FantasticFox1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

prop_38_noise <- noise_out %>% 
  mutate(Overall = ifelse(relerror_PCP_LOD < relerror_PCP_sqrt2, 1, 0),
         BelowLOD = ifelse(belowlod_relerror_PCP_LOD < belowlod_relerror_PCP_sqrt2, 1, 0),
         SafeRows = ifelse(safe_relerror_PCP_LOD < safe_relerror_PCP_sqrt2, 1, 0)) %>%
  select(rank, delta, iter, Overall, BelowLOD, SafeRows) %>% 
  pivot_longer(Overall:SafeRows) %>% 
  group_by(rank, delta, name) %>%
  summarize(n = n(),
            Better = sum(value)/n()) %>% ungroup() %>% 
  mutate(rank = as.factor(rank),
         sigma = "0.6",
         delta = as.numeric(delta))

prop_38_noise_10 <- noise_out %>% 
  mutate(relerror_PCP_LOD_10 = relerror_PCP_LOD - (0.1*relerror_PCP_LOD),
         belowlod_relerror_PCP_LOD_10 = belowlod_relerror_PCP_LOD - (0.1*belowlod_relerror_PCP_LOD),
         safe_relerror_PCP_LOD_10 = safe_relerror_PCP_LOD - (0.1*safe_relerror_PCP_LOD)) %>% 
  mutate(Overall = ifelse(relerror_PCP_LOD_10 < relerror_PCP_sqrt2, 1, 0),
         BelowLOD = ifelse(belowlod_relerror_PCP_LOD_10 < belowlod_relerror_PCP_sqrt2, 1, 0),
         SafeRows = ifelse(safe_relerror_PCP_LOD_10 < safe_relerror_PCP_sqrt2, 1, 0)) %>%
  select(rank, delta, iter, Overall, BelowLOD, SafeRows) %>% 
  pivot_longer(Overall:SafeRows) %>% 
  group_by(rank, delta, name) %>%
  summarize(n = n(),
            Better = sum(value)/n()) %>% ungroup() %>% 
  mutate(rank = as.factor(rank),
         sigma = "0.6",
         delta = as.numeric(delta))

prop_38_noise_5 <- noise_out %>% 
  mutate(relerror_PCP_LOD_10 = relerror_PCP_LOD - (0.05*relerror_PCP_LOD),
         belowlod_relerror_PCP_LOD_10 = belowlod_relerror_PCP_LOD - (0.05*belowlod_relerror_PCP_LOD),
         safe_relerror_PCP_LOD_10 = safe_relerror_PCP_LOD - (0.05*safe_relerror_PCP_LOD)) %>% 
  mutate(Overall = ifelse(relerror_PCP_LOD_10 < relerror_PCP_sqrt2, 1, 0),
         BelowLOD = ifelse(belowlod_relerror_PCP_LOD_10 < belowlod_relerror_PCP_sqrt2, 1, 0),
         SafeRows = ifelse(safe_relerror_PCP_LOD_10 < safe_relerror_PCP_sqrt2, 1, 0)) %>%
  select(rank, delta, iter, Overall, BelowLOD, SafeRows) %>% 
  pivot_longer(Overall:SafeRows) %>% 
  group_by(rank, delta, name) %>%
  summarize(n = n(),
            Better = sum(value)/n()) %>% ungroup() %>% 
  mutate(rank = as.factor(rank),
         sigma = "0.6",
         delta = as.numeric(delta))

prop_10 %>% ggplot(aes(y = Better, x = delta)) +
  geom_line(aes(group = rank, color = rank)) +
  geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed") +
  facet_wrap(~name) + scale_color_manual(values=wes_palette(n=5, name="FantasticFox1"))


noise_prop <- noise_out %>% 
  mutate(Overall = ifelse(relerror_PCP_LOD < relerror_PCP_sqrt2, 1, 0),
         BelowLOD = ifelse(belowlod_relerror_PCP_LOD < belowlod_relerror_PCP_sqrt2, 1, 0),
         SafeRows = ifelse(safe_relerror_PCP_LOD < safe_relerror_PCP_sqrt2, 1, 0)) %>%
  select(rank, delta, iter, Overall, BelowLOD, SafeRows) %>% 
  pivot_longer(Overall:SafeRows) %>% 
  group_by(rank, delta, name) %>%
  summarize(n = n(),
    Better = sum(value)/n()) %>% ungroup() %>% 
  mutate(rank = as.factor(rank),
         sigma = "0.6",
         delta = as.numeric(delta))

# Optimal MU
noise_prop %>% ggplot(aes(y = Better, x = delta)) +
  geom_line(aes(group = rank, color = rank)) +
  geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed") +
  facet_wrap(~name) + scale_color_manual(values=wes_palette(n=5, name="FantasticFox1"))


noise_plot %>% group_by(delta) %>%
  summarize(sd = sd(rel_error, na.rm = TRUE))
