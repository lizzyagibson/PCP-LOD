library(tidyverse)
library(gridExtra)
library(Matrix)
library(matconv)
library(patchwork)
library(janitor)
library(ggcorrplot)
library(ggfortify)  
library(factoextra)
library(knitr)
library(haven)
library(rlist)
library(mvtnorm)
library(reshape2)

# Load NHANES data
nhanes <- read_sas(here::here("./Data/studypop_lod.sas7bdat")) %>% clean_names()

prop <- function (x) {1 - (sum(x, na.rm = TRUE)/nrow(nhanes))}

# F08 <LOD
# Choose 10 chemicals with % detected > 80%
names <- nhanes %>% select(names(.)[grep("lc", names(.))]) %>% 
  summarize_all(prop) %>% select_if(~. > 0.8) %>% names() %>% str_sub(., 4, 6) %>% str_c("lbx", ., "la") %>% as.vector()

pops <- nhanes %>% 
  select(!!names) %>% na.omit(.)

## Rename
names(pops) <- str_sub(names(pops), 1, 6)
names(pops) <- str_replace(names(pops), "lbxd", "D")
names(pops) <- str_replace(names(pops), "lbxf", "F")
names(pops) <- str_replace(names(pops), "lbx", "PCB")
pops

## Vector of NHANES means
## log to approx normal dist
means <- as_vector(map(log(pops), function(x) mean(x, na.rm = TRUE)))

## Covariance matrix from NHANES
## log to approx normal dist
covs <- cov(log(pops))

## Simulate with multivariate normal function
## exp multi-normal to get multi-log normal
set.seed(1988)
sim <- exp(rmvnorm(100, mean = means, sigma = covs)) %>% as_tibble()
sim <- as_tibble(apply(sim, 2, scale, center = FALSE))

## Viz
sim %>%
  mutate(id = 1:nrow(.)) %>% 
  gather(key = pop, value = value, -id) %>% 
  ggplot(aes(x = value)) +
  geom_histogram() + facet_wrap(~pop) +
  theme_minimal()

# write_csv(sim, "mix_data_lod_0.csv")

## Push <LOD to -1

### Create version with 10% lowest values for F08 as below the LOD
mix_data_lod_10 <- sim %>% mutate(F08 = ifelse(F08 < quantile(F08, probs = .10), -1, F08))
# write_csv(mix_data_lod_10, "mix_data_lod_10.csv")

### Create version with 20% lowest values for F08 as below the LOD
mix_data_lod_20 <- sim %>% mutate(F08 = ifelse(F08 < quantile(F08, probs = .20), -1, F08)) 
# write_csv(mix_data_lod_20, "mix_data_lod_20.csv")

### Create version with 30% lowest values for F08 as below the LOD
mix_data_lod_30 <- sim %>% mutate(F08 = ifelse(F08 < quantile(F08, probs = .30), -1, F08)) 
# write_csv(mix_data_lod_30, "mix_data_lod_30.csv")

# Create version with 40% lowest values for F08 as below the LOD
mix_data_lod_40 <- sim %>% mutate(F08 = ifelse(F08 < quantile(F08, probs = .40), -1, F08)) 
# write_csv(mix_data_lod_40, "mix_data_lod_40.csv")

### Create version with 50% lowest values for F08 as below the LOD
mix_data_lod_50 <- sim %>% mutate(F08 = ifelse(F08 < quantile(F08, probs = .50), -1, F08)) 
# write_csv(mix_data_lod_50, "mix_data_lod_50.csv")

### Create quantile LOD lists
delta10 <- c(rep(0, times = 6), quantile(sim$F08, probs = 0.10), 0, 0)
delta20 <- c(rep(0, times = 6), quantile(sim$F08, probs = 0.20), 0, 0)
delta30 <- c(rep(0, times = 6), quantile(sim$F08, probs = 0.30), 0, 0)
delta40 <- c(rep(0, times = 6), quantile(sim$F08, probs = 0.40), 0, 0)
delta50 <- c(rep(0, times = 6), quantile(sim$F08, probs = 0.50), 0, 0)

# All chemicals <LOD
# Choose chemicals detected > 60%
names_all <- nhanes %>% select(names(.)[grep("lc", names(.))]) %>% 
  summarize_all(prop) %>% select_if(~. > 0.6) %>% names() %>% str_sub(., 4, 6) %>% str_c("lbx", ., "la") %>% as.vector()

## Select chemicals
pops_all <- nhanes %>% select(!!names_all) %>% na.omit(.)

## Rename
names(pops_all) <- str_sub(names(pops_all), 1, 6)
names(pops_all) <- str_replace(names(pops_all), "lbxd", "D")
names(pops_all) <- str_replace(names(pops_all), "lbxf", "F")
names(pops_all) <- str_replace(names(pops_all), "lbx", "PCB")
pops_all

## Vector of NHANES means
## log to approx normal dist
means_all <- as_vector(map(log(pops_all), function(x) mean(x, na.rm = TRUE)))

## Covariance matrix from NHANES
## log to approx normal dist
covs_all <- cov(log(pops_all))

## Simulate with multivariate normal function
## exp multi-normal to get multi-log normal
set.seed(1988)
sim_all <- exp(rmvnorm(1000, mean = means_all, sigma = covs_all)) %>% as_tibble()

## Scale Simulations
## Divide by standard deviation, do not mean center.  
sim_all <- as_tibble(apply(sim_all, 2, scale, center = FALSE))


## Viz
sim_all %>%
  mutate(id = 1:nrow(.)) %>% 
  gather(key = pop, value = value, -id) %>% 
  ggplot(aes(x = value)) +
  geom_histogram() + facet_wrap(~pop) +
  theme_minimal()

# write_csv(sim_all, "mix_data_lod_0_all.csv")

## Create \<LOD Datasets

## Push to -1

### Create version with 10% lowest values for each variable as below the LOD
mix_data_lod_10_all <- sim_all %>% mutate_all(~ifelse(. < quantile(., probs = .10), -1, .)) 
# write_csv(mix_data_lod_10_all, "mix_data_lod_10_all.csv")

### Create version with 20% lowest values for each variable as below the LOD
mix_data_lod_20_all <- sim_all %>%  mutate_all(~ifelse(. < quantile(., probs = .20), -1, .)) 
# write_csv(mix_data_lod_20_all, "mix_data_lod_20_all.csv")

### Create version with 30% lowest values for each variable as below the LOD
mix_data_lod_30_all <- sim_all %>% mutate_all(~ifelse(. < quantile(., probs = .30), -1, .)) 
# write_csv(mix_data_lod_30_all, "mix_data_lod_30_all.csv")

### Create version with 40% lowest values for each variable as below the LOD
mix_data_lod_40_all <- sim_all %>% mutate_all(~ifelse(. < quantile(., probs = .40), -1, .)) 
# write_csv(mix_data_lod_40_all, "mix_data_lod_40_all.csv")

### Create version with 50% lowest values for each variable as below the LOD
mix_data_lod_50_all <- sim_all %>%  mutate_all(~ifelse(. < quantile(., probs = .50), -1, .)) 
# write_csv(mix_data_lod_50_all, "mix_data_lod_50_all.csv")

## Quantiles = LOD
delta10_all <- sim_all %>% summarise_all(quantile, probs = .10) %>% as_vector()
delta20_all <- sim_all %>% summarise_all(quantile, probs = .20) %>% as_vector()
delta30_all <- sim_all %>% summarise_all(quantile, probs = .30) %>% as_vector()
delta40_all <- sim_all %>% summarise_all(quantile, probs = .40) %>% as_vector()
delta50_all <- sim_all %>% summarise_all(quantile, probs = .50) %>% as_vector()
