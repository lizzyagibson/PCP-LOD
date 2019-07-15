library(R.matlab)
library(tidyverse)

mixture <- readMat("./MATLAB/mixtures_data.mat")

mixture_data <- as.data.frame(mixture) %>% as_tibble() %>% 
  select(Al, As, Ba, bc, Br, Ca, Cl,
         Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
         Ti,  V, Zn) %>% 
  drop_na() %>% as.matrix()

# Create version with 10% below the LOD
mix_data_lod_10 <- mixture_data

set.seed(8)
mysample <- sample(length(mix_data_lod_10), size=0.1*ncol(mix_data_lod_10)*nrow(mix_data_lod_10), replace =F)
length(mysample)

mix_data_lod_10[mysample] <- -1

summary(mix_data_lod_10)
mix_data_lod_10[mix_data_lod_10 == -1]
head(mix_data_lod_10)

write_csv(as_tibble(mix_data_lod_10), "./Below_LOD/mix_data_lod_10.csv")

