library(R.matlab)
library(tidyverse)
source("./Original/R Code/stable_pcp_alternating.R")

mixture <- readMat("./Data/mixtures_data.mat")

L_matlab <- read_csv("./Original/MATLAB/L_matlab.csv", col_names = FALSE)
S_matlab <- read_csv("./Original/MATLAB/S_matlab.csv", col_names = FALSE) %>% as.matrix()

mixture_data <- as.data.frame(mixture) %>% as_tibble() %>% select(Al, As, Ba, bc, Br, Ca, Cl,
                                                                  Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
                                                                  Ti,  V, Zn)

mix <- mixture_data[complete.cases(mixture_data),]

X <- scale(mix, center = TRUE, scale = TRUE)
#should we always normalize data beforehand?

m <- nrow(X)
n <- ncol(X)

lambda = 1/sqrt(m)

mixture_out <- stable_pcp_alternating(X, 4/sqrt(m), 10)

mixture_S <- mixture_out$S
mixture_L <- mixture_out$L

svd(mixture_L)$d
svd(X)$d

#Compare with MATLAB output
colSums(mixture_S)
colSums(S_matlab)

svd(mixture_L)$d

colSums(mixture_L)
colSums(L_matlab)

#compare svd(L) with MATLAB code
U_matlab <- read_csv("./U_matlab.csv", col_names = FALSE)
sigma_matlab <- read_csv("./Sigma_matlab.csv", col_names = FALSE) %>% as.matrix()
V_matlab <- read_csv("./V_matlab.csv", col_names = FALSE)

#similar sigma!
colSums(sigma_matlab)
svd(mixture_L)$d

#similar V, too!
colSums(V_matlab)
colSums(svd(mixture_L)$v)

summary(V_matlab$X1)
summary(svd(mixture_L)$v[,1])

hist(V_matlab$X1)
hist(svd(mixture_L)$v[,1])
