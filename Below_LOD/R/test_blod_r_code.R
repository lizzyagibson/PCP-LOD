source("/Users/lizzy/Principle.Component.Pursuit/Below_LOD/R/function_blod.R")

library(R.matlab)
library(tidyverse)

mixture <- readMat("./Data/mixtures_data.mat")

mixture_data <- as.data.frame(mixture) %>% as_tibble() %>% select(Al, As, Ba, bc, Br, Ca, Cl,
                                                                  Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
                                                                  Ti,  V, Zn)

mix <- mixture_data[complete.cases(mixture_data),]

X <- scale(mix, center = TRUE, scale = TRUE)

m <- nrow(X)
n <- ncol(X)

lambda = 1/sqrt(m)

mixture_out <- pcp_lod(X, 4/sqrt(m), 10, 0)
