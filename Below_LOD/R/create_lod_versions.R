library(R.matlab)
library(tidyverse)

# Read air pollution data
mixture <- readMat("./Data/mixtures_data.mat")

mixture_data <- as.data.frame(mixture) %>% as_tibble() %>% 
  select(Al, As, Ba, bc, Br, Ca, Cl,
         Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
         Ti,  V, Zn) %>% 
  drop_na()

# Create version with 10% lowest values for each variable as below the LOD
ifelse(mixture_data$Al <= quantile(mixture_data$Al, probs = .10), -1, mixture_data$Al)
sum(mixture_data$Al <= quantile(mixture_data$Al, probs = .10))/nrow(mixture_data)

mix_data_lod_10 <- mixture_data %>% 
  mutate_all(~ifelse(. <= quantile(., probs = .10), -1, .))

summary(mix_data_lod_10)
#write_csv(as_tibble(mix_data_lod_10), "./Below_LOD/R/BLOD_airpol_data/mix_data_lod_10.csv")

# Create version with 20% lowest values for each variable as below the LOD
ifelse(mixture_data$Al <= quantile(mixture_data$Al, probs = .20), -1, mixture_data$Al)
sum(mixture_data$Al <= quantile(mixture_data$Al, probs = .20))/nrow(mixture_data)

mix_data_lod_20 <- mixture_data %>% 
  mutate_all(~ifelse(. <= quantile(., probs = .20), -1, .))

summary(mix_data_lod_20)
#write_csv(as_tibble(mix_data_lod_20), "./Below_LOD/R/BLOD_airpol_data/mix_data_lod_20.csv")

# Create version with 30% lowest values for each variable as below the LOD
ifelse(mixture_data$Al <= quantile(mixture_data$Al, probs = .30), -1, mixture_data$Al)
sum(mixture_data$Al <= quantile(mixture_data$Al, probs = .30))/nrow(mixture_data)

mix_data_lod_30 <- mixture_data %>% 
  mutate_all(~ifelse(. <= quantile(., probs = .30), -1, .))

summary(mix_data_lod_30)
#write_csv(as_tibble(mix_data_lod_30), "./Below_LOD/R/BLOD_airpol_data/mix_data_lod_30.csv")

# Create version with 40% lowest values for each variable as below the LOD
ifelse(mixture_data$Al <= quantile(mixture_data$Al, probs = .40), -1, mixture_data$Al)
sum(mixture_data$Al <= quantile(mixture_data$Al, probs = .40))/nrow(mixture_data)

mix_data_lod_40 <- mixture_data %>% 
  mutate_all(~ifelse(. <= quantile(., probs = .40), -1, .))

summary(mix_data_lod_40)
#write_csv(as_tibble(mix_data_lod_40), "./Below_LOD/R/BLOD_airpol_data/mix_data_lod_40.csv")

# Create version with 50% lowest values for each variable as below the LOD
ifelse(mixture_data$Al <= quantile(mixture_data$Al, probs = .50), -1, mixture_data$Al)
sum(mixture_data$Al <= quantile(mixture_data$Al, probs = .50))/nrow(mixture_data)

mix_data_lod_50 <- mixture_data %>% 
  mutate_all(~ifelse(. <= quantile(., probs = .50), -1, .))

summary(mix_data_lod_50)
#write_csv(as_tibble(mix_data_lod_50), "./Below_LOD/R/BLOD_airpol_data/mix_data_lod_50.csv")


## Quantiles

mixture_data %>% 
  summarise_all(quantile, probs = .10) %>% as.matrix()

mixture_data %>% 
  summarise_all(quantile, probs = .20) %>% as.matrix()

mixture_data %>% 
  summarise_all(quantile, probs = .30) %>% as.matrix()

mixture_data %>% 
  summarise_all(quantile, probs = .40) %>% as.matrix()

mixture_data %>% 
  summarise_all(quantile, probs = .50) %>% as.matrix()
