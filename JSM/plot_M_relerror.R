# Plot PCP-LOD, LOD/sqrt(2), PCA resulting error
# Across rank, sigma, mu, and LOD
library(tidyverse)
library(R.matlab)
library(RColorBrewer)
library(wesanderson)
theme_set(theme_minimal() + theme(legend.position = "bottom"))

rank = c(1, 2)
sigma = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, "1.0")
delta = c(0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

# r = 1
# d = 0.3
# s = 0.1

pcp_out <- tibble()

for (r in rank) {
for (d in delta) {
for (s in sigma) {
#for (r in rank) {
  # Mu without sigma
  
  read <- readMat(here::here(paste0("HPC_PCP/Var/data_lod_", d, "_rank_", r, 
                                    "_sigma_", s, ".mat")))
  
  relerror <- read[[1]] %>% 
    as_tibble() %>% rename("relerror_PCA" =1, "relerror_PCP_sqrt2" =2, "relerror_PCP_LOD" =3)
  
  # score_error <- read[[3]] %>% 
  #   as_tibble() %>% rename("score_relerror_PCA" =1, "score_relerror_PCP_sqrt2" =2, 
  #                          "score_relerror_PCP_LOD" =3) 
  # 
  # pattern_error <- read[[4]] %>% 
  #   as_tibble() %>% rename("pattern_relerror_PCA" =1, "pattern_relerror_PCP_sqrt2" =2, 
  #                          "pattern_relerror_PCP_LOD" =3) 
  out_mu <- read[[2]] %>% 
    as_tibble() %>% rename("belowlod_relerror_PCA" =1, "belowlod_relerror_PCP_sqrt2" =2, 
                           "belowlod_relerror_PCP_LOD" =3) %>% 
    mutate(rank = r,
           sigma = s,
           delta = d,
           iter = 1:100) %>% 
    cbind(., relerror) %>% as_tibble()
  
  
  pcp_out <- rbind(pcp_out, out_mu) %>% drop_na(.)
  }
}
}

pcp_plot <- pcp_out %>% 
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
  select(-name)

names(pcp_plot)

options(
  ggplot2.continuous.colour = "dark2",
  ggplot2.continuous.fill = "dark2"
)

scale_colour_discrete = scale_colour_viridis_d
scale_fill_discrete = scale_fill_viridis_d

# rank increase -- PCP_LOD does better
pcp_plot %>% 
  filter(rank == "1") %>% 
  filter(model != "PCA") %>% 
  #filter(type != "overall") %>% 
  #filter(type != "belowlod") %>% 
  filter(sigma == "0.1") %>% 
  ggplot(aes(x = delta, y = rel_error, color = model)) +
  geom_boxplot(notch = FALSE, outlier.shape = NA) +
  facet_grid(~type) + scale_y_log10() +
  theme_bw() + scale_color_manual(values=wes_palette(n=3, name="FantasticFox1"))
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")





