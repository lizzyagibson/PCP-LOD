library(R.matlab)
library(tidyverse)

mixture <- readMat("/Users/lizzy/Principle.Component.Pursuit/Data/mixtures_data.mat")

mixture_data <- as.data.frame(mixture) %>% 
  as_tibble() %>% dplyr::select(Al, As, Ba, bc, Br, Ca, Cl,
                         Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
                         Ti,  V, Zn) %>% as.matrix()

X <- mixture_data[complete.cases(mixture_data),]
X[1:5, 1:5]

m <- nrow(X)
n <- ncol(X)

root_out = root_pcp(X, 1/sqrt(m), 10)

norm(root_out[[1]], "F")
norm(root_out[[2]], "F")
