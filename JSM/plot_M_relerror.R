# Plot PCP-LOD, LOD/sqrt(2), PCA resulting error
# Across rank, sigma, mu, and LOD
library(tidyverse)
library(R.matlab)
library(RColorBrewer)
library(wesanderson)
theme_set(theme_bw(base_size = 20) + theme(legend.position = "bottom",
                             strip.background =element_rect(fill="white")))

rank = c(1, 2, 3, 4, 5)
sigma = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
delta = c(0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

dim(expand.grid(rank, sigma, delta))

# r = 1
# d = 0.3
# s = 0.1

pcp_out <- tibble()

for (r in rank) {
for (d in delta) {
for (s in sigma) {
  # Mu without sigma
  
  read <- readMat(here::here(paste0("HPC_PCP/pcp_all/pcp_all_lod_", d, "_rank_", r, 
                                    "_sigma_", s, ".mat")))
  
  relerror <- read[[1]] %>% 
    as_tibble() %>% rename("relerror_PCA" =1, "relerror_PCP_sqrt2" =2, "relerror_PCP_LOD" =3)
  
  score_error <- read[[3]] %>%
    as_tibble() %>% rename("score_relerror_PCA" =1, "score_relerror_PCP_sqrt2" =2,
                           "score_relerror_PCP_LOD" =3)

  pattern_error <- read[[4]] %>%
    as_tibble() %>% rename("pattern_relerror_PCA" =1, "pattern_relerror_PCP_sqrt2" =2,
                           "pattern_relerror_PCP_LOD" =3)
  out_mu <- read[[2]] %>% 
    as_tibble() %>% rename("belowlod_relerror_PCA" =1, "belowlod_relerror_PCP_sqrt2" =2, 
                           "belowlod_relerror_PCP_LOD" =3) %>% 
    mutate(rank = r,
           sigma = s,
           delta = d,
           iter = 1:100) %>% 
    cbind(., relerror, score_error, pattern_error) %>% as_tibble()
  
  
  pcp_out <- rbind(pcp_out, out_mu) %>% drop_na(.)
  }
}
}

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


names(pcp_plot)

# \definecolor{matblue}{HTML}{2196F3}
# \definecolor{matbluelight}{HTML}{BBDEFB}
# \definecolor{matbluedark}{HTML}{1976D2}

# rank increase -- PCP_LOD does better
jsm_overall <- 
  pcp_plot %>% 
  filter(rank == "4") %>% 
  filter(sigma == "0.6") %>% 
  filter(delta != "1") %>% 
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

# FantasticFox1 = c("#E58601", "#B40F20"),  
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("/Users/lizzy/Summer 2020/JSM PCP-LOD/figures/jsm_overall.pdf", width = 11)
jsm_overall
dev.off()

##
# Loadings and scores
##

svd_plot <- pcp_out %>%
  select(rank, sigma, delta, iter, score_relerror_PCA:score_relerror_PCP_LOD, 
         pattern_relerror_PCA:pattern_relerror_PCP_LOD) %>% 
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
  mutate(type = case_when(type == "score" ~ "Scores",
                          type == "pattern" ~ "Loadings"),
         type = fct_relevel(type, "Overall", "< LOD")) %>% 
  mutate(model = case_when(model == "PCA" ~ "PCA",
                           model == "PCP_LOD" ~ "PCP-LOD",
                           model == "PCP_sqrt2" ~ #paste0("PCP w/ LOD/", intToUtf8(8730),2))) %>% 
                             "PCP"))


names(pcp_plot)

# rank increase -- PCP_LOD does better
jsm_svd <- 
svd_plot %>% 
  filter(rank == "4") %>% 
  filter(sigma == "0.6") %>% 
  filter(delta != "1") %>% 
  ggplot(aes(x = delta, y = rel_error, color = model)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.25) +
  facet_wrap(~type, scales = "free_y") + scale_y_log10() +
  scale_color_manual(values = c("#E58601", # PCA
                                "#B40F20",  # PCP SQRT2
                                "#2196F3" # PCP LOD
  )) +
  labs(x = "Values < LOD (proportion)",
       y = "Relative Prediction Error") +
  theme(
    legend.title = element_blank())

# FantasticFox1 = c("#E58601", "#B40F20"),  
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("/Users/lizzy/Summer 2020/JSM PCP-LOD/figures/jsm_svd.pdf", width = 11)
jsm_svd
dev.off()


