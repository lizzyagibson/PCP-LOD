library(pcpr)
library(tidyverse)
library(janitor)
library(tictoc)

metab <- read_csv("Data/metabolomics.data.csv") %>% clean_names() %>% select(-x1) %>% as.matrix()
head(metab)[,1:10]
dim(metab)
summary(metab[,1:6])
apply(metab[,1:10], 2, sd)

# flip long and short dims
mm <- nrow(metab)
nn <- ncol(metab)

mu <- sqrt(mm/(2*log(nn*mm)))
lam <- 1/sqrt(nn)

tic()
raw_out <- pcp_lod(metab, lam,  mu , 0)
toc()
# Did not converge
# 2.5 hours

# log 2
logged <- as_tibble(metab) %>% 
  mutate_all(function (x) ifelse(x == 0, x + 1, x)) %>% 
  mutate_all(log2) %>% as.matrix()

tic()
log_out <- pcp_lod(logged, lam,  mu , 0)
toc()
# 88.122 sec elapsed

# std
stnd <- as_tibble(metab) %>% 
  mutate_all(function(x) x/sd(x)) %>% 
  as.matrix()

tic()
std_out <- pcp_lod(stnd, lam,  mu , 0)
toc()
# 47.952 sec elapsed
