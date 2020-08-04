# Plot PCP-LOD, LOD/sqrt(2), PCA resulting error
# Across rank, sigma, mu, and LOD
library(tidyverse)
library(R.matlab)
library(RColorBrewer)
library(wesanderson)
theme_set(theme_bw(base_size = 20) + theme(legend.position = "bottom",
                             strip.background =element_rect(fill="white")))

rank = c(1, 2, 3, 4, 5)
delta = c("0", 0.05, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

# r = 1
# d = 0.3

nonoise_38 <- tibble()

for (r in rank) {
  for (d in delta) {
    read <- readMat(here::here(paste0("HPC_PCP/nonoise/no_noise_38_lod_", d, "_rank_", r, ".mat")))
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
    nonoise_38 <- rbind(nonoise_38, out) %>% drop_na(iter)
  }
}

plot_38 <- 
  nonoise_38 %>% 
  pivot_longer(grep("relerror", colnames(.)),
               values_to = "rel_error") %>% 
  mutate(type = ifelse(grepl("belowlod", name), "belowlod", 
                       ifelse(grepl("safe", name), "safe", "overall")),
         model = ifelse(grepl("PCA", name), "PCA",
                        ifelse(grepl("sqrt2", name), "PCP_sqrt2", "PCP_LOD")),
         rank = as.factor(rank)) %>% 
  select(-name)

prop_38 <- nonoise_38 %>% 
  mutate(Overall = ifelse(relerror_PCP_LOD < relerror_PCP_sqrt2, 1, 0),
         BelowLOD = ifelse(belowlod_relerror_PCP_LOD < belowlod_relerror_PCP_sqrt2, 1, 0),
         SafeRows = ifelse(safe_relerror_PCP_LOD < safe_relerror_PCP_sqrt2, 1, 0)) %>%
  select(rank, delta, iter, Overall, BelowLOD, SafeRows) %>% 
  pivot_longer(Overall:SafeRows) %>% 
  group_by(rank, delta, name) %>%
  summarize(n = n(),
            Better = sum(value)/n()) %>% ungroup() %>% 
  mutate(rank = as.factor(rank))

prop_38_5 <- nonoise_38 %>% 
  mutate(relerror_PCP_LOD_10 = relerror_PCP_LOD - (0.05*relerror_PCP_LOD),
        belowlod_relerror_PCP_LOD_10 = belowlod_relerror_PCP_LOD - (0.05*belowlod_relerror_PCP_LOD),
        safe_relerror_PCP_LOD_10 = safe_relerror_PCP_LOD - (0.05*safe_relerror_PCP_LOD)) %>% 
  mutate(Overall = ifelse(relerror_PCP_LOD_10 < relerror_PCP_sqrt2, 1, 0),
         BelowLOD = ifelse(belowlod_relerror_PCP_LOD_10 < belowlod_relerror_PCP_sqrt2, 1, 0),
         SafeRows = ifelse(safe_relerror_PCP_LOD_10 < safe_relerror_PCP_sqrt2, 1, 0)) %>%
  # mutate(Overall = ifelse(relerror_PCP_LOD < relerror_PCP_sqrt2, 1, 0),
  #        BelowLOD = ifelse(belowlod_relerror_PCP_LOD < belowlod_relerror_PCP_sqrt2, 1, 0),
  #        SafeRows = ifelse(safe_relerror_PCP_LOD < safe_relerror_PCP_sqrt2, 1, 0)) %>%
  select(rank, delta, iter, Overall, BelowLOD, SafeRows) %>% 
  pivot_longer(Overall:SafeRows) %>% 
  group_by(rank, delta, name) %>%
  summarize(n = n(),
    Better = sum(value)/n()) %>% ungroup() %>% 
  mutate(rank = as.factor(rank))

prop_38_10 <- nonoise_38 %>% 
  mutate(relerror_PCP_LOD_10 = relerror_PCP_LOD - (0.1*relerror_PCP_LOD),
         belowlod_relerror_PCP_LOD_10 = belowlod_relerror_PCP_LOD - (0.1*belowlod_relerror_PCP_LOD),
         safe_relerror_PCP_LOD_10 = safe_relerror_PCP_LOD - (0.1*safe_relerror_PCP_LOD)) %>% 
  mutate(Overall = ifelse(relerror_PCP_LOD_10 < relerror_PCP_sqrt2, 1, 0),
         BelowLOD = ifelse(belowlod_relerror_PCP_LOD_10 < belowlod_relerror_PCP_sqrt2, 1, 0),
         SafeRows = ifelse(safe_relerror_PCP_LOD_10 < safe_relerror_PCP_sqrt2, 1, 0)) %>%
  # mutate(Overall = ifelse(relerror_PCP_LOD < relerror_PCP_sqrt2, 1, 0),
  #        BelowLOD = ifelse(belowlod_relerror_PCP_LOD < belowlod_relerror_PCP_sqrt2, 1, 0),
  #        SafeRows = ifelse(safe_relerror_PCP_LOD < safe_relerror_PCP_sqrt2, 1, 0)) %>%
  select(rank, delta, iter, Overall, BelowLOD, SafeRows) %>% 
  pivot_longer(Overall:SafeRows) %>% 
  group_by(rank, delta, name) %>%
  summarize(n = n(),
            Better = sum(value)/n()) %>% ungroup() %>% 
  mutate(rank = as.factor(rank))

#jsm_prop <- 
prop_38 %>%
  #filter(delta != "0") %>% 
  mutate(name = case_when(name == "BelowLOD" ~ "< LOD",
                          name == "Overall" ~ "Overall",
                          name == "SafeRows" ~ "Safe Rows"),
         name = fct_relevel(name, "Overall", "< LOD", "Safe Rows")) %>% 
  mutate(sigma = as.factor("0.0")) %>% 
  filter(rank == "4") %>% 
  ggplot(aes(y = Better, x = delta, color = sigma)) +
  geom_line(aes(group = 1), size = 1.1) +
  geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed") +
  facet_wrap(~name) + 
  scale_y_continuous(breaks=seq(0,1.1,0.25)) +
  labs(x = "Values < LOD (proportion)",
       y = "Better performance by PCP-LOD (proportion)",
       color = "Noise level:") +
  scale_color_manual(values = c("#2196F3"))

#pdf("/Users/lizzy/Summer 2020/JSM PCP-LOD/figures/jsm_prop.pdf", width = 12)
jsm_prop
#dev.off()

compare <- prop_38 %>%
  mutate(sigma = as.factor("0.0")) %>% 
  rbind(., prop_38_noise) %>% 
  #filter(delta != "0") %>% 
  mutate(name = case_when(name == "BelowLOD" ~ "< LOD",
                          name == "Overall" ~ "Overall",
                          name == "SafeRows" ~ "Safe Rows"),
         name = fct_relevel(name, "Overall", "< LOD", "Safe Rows")) %>% 
  filter(rank == "4") %>% 
  ggplot(aes(y = Better, x = delta, color = sigma)) +
  geom_line(aes(group = sigma), size = 1.1) +
  geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed") +
  facet_wrap(~name) + 
  scale_y_continuous(breaks=seq(0,1.1,0.25)) +
  labs(x = "Values < LOD (proportion)",
       y = "Better performance by PCP-LOD (proportion)",
       color = "Noise level:") +
  scale_color_manual(values = c("#2196F3",
                                "#FD6467"))
                                #"#81A88D")) 

compare_10 <- prop_38_10 %>%
  mutate(sigma = as.factor("0.0")) %>% 
  rbind(., prop_38_noise_10) %>% 
  #filter(delta != "0") %>% 
  mutate(name = case_when(name == "BelowLOD" ~ "< LOD",
                          name == "Overall" ~ "Overall",
                          name == "SafeRows" ~ "Safe Rows"),
         name = fct_relevel(name, "Overall", "< LOD", "Safe Rows")) %>% 
  filter(rank == "4") %>% 
  ggplot(aes(y = Better, x = delta, color = sigma)) +
  geom_line(aes(group = sigma), size = 1.1) +
  geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed") +
  facet_wrap(~name) + 
  scale_y_continuous(breaks=seq(0,1.1,0.25)) +
  labs(x = "Values < LOD (proportion)",
       y = "Better performance by PCP-LOD (proportion)",
       color = "Noise level:") +
  scale_color_manual(values = c("#2196F3",
                                "#FD6467"))

compare_5 <- prop_38_5 %>%
  mutate(sigma = as.factor("0.0")) %>% 
  rbind(., prop_38_noise_5) %>% 
  #filter(delta != "0") %>% 
  mutate(name = case_when(name == "BelowLOD" ~ "< LOD",
                          name == "Overall" ~ "Overall",
                          name == "SafeRows" ~ "Safe Rows"),
         name = fct_relevel(name, "Overall", "< LOD", "Safe Rows")) %>% 
  filter(rank == "4") %>% 
  ggplot(aes(y = Better, x = delta, color = sigma)) +
  geom_line(aes(group = sigma), size = 1.1) +
  geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed") +
  facet_wrap(~name) + 
  scale_y_continuous(breaks=seq(0,1.1,0.25)) +
  labs(x = "Values < LOD (proportion)",
       y = "Better performance by PCP-LOD (proportion)",
       color = "Noise level:") +
  scale_color_manual(values = c("#2196F3",
                                "#FD6467"))

pdf("/Users/lizzy/Summer 2020/JSM PCP-LOD/figures/compare.pdf", width = 12)
compare
dev.off()

pdf("/Users/lizzy/Summer 2020/JSM PCP-LOD/figures/compare10.pdf", width = 12)
compare_10
dev.off()

pdf("/Users/lizzy/Summer 2020/JSM PCP-LOD/figures/compare5.pdf", width = 12)
compare_5
dev.off()
