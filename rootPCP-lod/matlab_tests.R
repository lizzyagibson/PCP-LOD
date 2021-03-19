# File: matlab_tests.R
# Author: Lawrence Chillrud (lgc2139@cumc.columbia.edu)
# Date: 3/16/21

# This file details non-convex PCP tests we would like to attempt in matlab,
# to rule out transcription error in the R versions of non-convex PCP. 

######## 0. PACKAGE IMPORTS ######## 

# for PCP
library(PCPhelpers)
library(pcpr)

# for data wrangling
library(tidyverse)

# for data viz.
library(heatmaply)
library(ggplot2)
library(kableExtra)

######## 1. READ IN NYC DATA ######## 
prep_data <- function(
  file = here::here("./experiments/Data/nyc_datasets/nyc_daily_pm2.5_components.csv"), 
  years = c(as.Date("2007-01-01"), as.Date("2016-01-01")),
  exclude = c("cobalt", "antimony"),
  title = "NYC") 
{
  data <- readr::read_csv(file) %>% 
    filter(Date >= years[1], Date < years[2]) %>% 
    select(!all_of(exclude))
  
  data.mat <- data %>% select(!Date)
  data.mat.scaled <- scale(data.mat, center = FALSE)
  
  title <- paste0(title, ": ", format.Date(years[1], "%Y"), " - ", format.Date(years[2], "%Y"))
  
  list(M = data.mat, M.scaled = data.mat.scaled, M.dates = data$Date, M.title = title)
}

nyc <- prep_data()

######## 2. ASSIGN MAT OF INTEREST ######## 
mat <- nyc$M.scaled %>% as_tibble() %>% mutate(Date = nyc$M.dates) %>% na.omit()
dates <- mat$Date
mat <- mat %>% select(!Date) %>% as.matrix() 

summary <- data.frame() # for summary statistics later...

######## 3. CORRUPT MAT ######## 
cor_mat <- corrupt_mat_randomly(mat, k = 1, perc_b = 0.2) # a list containing cor.mat and cor.mask
mask <- cor_mat$cor.mask

# send to matlab
# write_csv(as.data.frame(mask), file = "./root_pcp/law_mask.csv")
# write_csv(as.data.frame(mat), file = "./root_pcp/law_mat.csv")
# write_csv(as.data.frame(cor_mat$cor.mat), file = "./root_pcp/law_cor_mat.csv")

######## 4. NONCVX LOD - GOOD RUN ######## 
lambda <- 0.122
mu <- 17.6
rank <- 6

pcp_out_lod_good <- root_pcp_noncvx_nonnegL_na_lod(cor_mat$cor.mat, lambda = lambda, mu = mu, r = rank, LOD = 0, verbose = T)

score <- norm((mat - pcp_out_lod_good$L - pcp_out_lod_good$S)*mask, "F") / norm(mat * mask, "F")
L.rank <- Matrix::rankMatrix(pcp_out_lod_good$L, tol = 1e-04)
S.sparsity <- sparsity(pcp_out_lod_good$S, tol = 1e-04)

Lnorm = norm(pcp_out_lod_good$L, "F")
Snorm = norm(pcp_out_lod_good$S, "F")

summary <- rbind(summary, data.frame(Algo = "NONCVX LOD - GOOD RUN", lambda = lambda, mu = mu, 
                                     rank = rank, Rel.Err = score, L.rank = L.rank, S.sparsity = S.sparsity,
                                     Lnorm = Lnorm, Snorm = Snorm))

######## 4. NONCVX LOD - BAD RUN ######## 
lambda <- 0.272
mu <- 1.61
rank <- 9

pcp_out_lod_bad <- root_pcp_noncvx_nonnegL_na_lod(cor_mat$cor.mat, lambda = lambda, mu = mu, r = rank, LOD = 0, verbose = T)

score <- norm((mat - pcp_out_lod_bad$L - pcp_out_lod_bad$S)*mask, "F") / norm(mat * mask, "F")
L.rank <- Matrix::rankMatrix(pcp_out_lod_bad$L, tol = 1e-04)
S.sparsity <- sparsity(pcp_out_lod_bad$S, tol = 1e-04)

Lnorm = norm(pcp_out_lod_bad$L, "F")
Snorm = norm(pcp_out_lod_bad$S, "F")

summary <- rbind(summary, data.frame(Algo = "NONCVX LOD - BAD RUN", lambda = lambda, mu = mu, 
                                     rank = rank, Rel.Err = score, L.rank = L.rank, S.sparsity = S.sparsity,
                                     Lnorm = Lnorm, Snorm = Snorm))

######## 5. REG NONCVX - GOOD RUN ######## 
lambda <- 0.272
mu <- 13.6
rank <- 4

pcp_out_good <- root_pcp_noncvx_nonnegL_na(cor_mat$cor.mat, lambda = lambda, mu = mu, r = rank, verbose = T)

score <- norm((mat - pcp_out_good$L - pcp_out_good$S)*mask, "F") / norm(mat * mask, "F")
L.rank <- Matrix::rankMatrix(pcp_out_good$L, tol = 1e-04)
S.sparsity <- sparsity(pcp_out_good$S, tol = 1e-04)

Lnorm = norm(pcp_out_good$L, "F")
Snorm = norm(pcp_out_good$S, "F")

summary <- rbind(summary, data.frame(Algo = "REG NONCVX - GOOD RUN", lambda = lambda, mu = mu, 
                                     rank = rank, Rel.Err = score, L.rank = L.rank, S.sparsity = S.sparsity,
                                     Lnorm = Lnorm, Snorm = Snorm))

######## 6. REG NONCVX - BAD RUN ######## 
lambda <- 0.182
mu <- 9.61
rank <- 6

pcp_out_bad <- root_pcp_noncvx_nonnegL_na(cor_mat$cor.mat, lambda = lambda, mu = mu, r = rank, verbose = T)

score <- norm((mat - pcp_out_bad$L - pcp_out_bad$S)*mask, "F") / norm(mat * mask, "F")
L.rank <- Matrix::rankMatrix(pcp_out_bad$L, tol = 1e-04)
S.sparsity <- sparsity(pcp_out_bad$S, tol = 1e-04)

Lnorm = norm(pcp_out_bad$L, "F")
Snorm = norm(pcp_out_bad$S, "F")

summary <- rbind(summary, data.frame(Algo = "REG NONCVX - BAD RUN", lambda = lambda, mu = mu, 
                                     rank = rank, Rel.Err = score, L.rank = L.rank, S.sparsity = S.sparsity,
                                     Lnorm = Lnorm, Snorm = Snorm))

######## 7. SUMMARY OF RESULTS ######## 
# see output on console for convergence per run; should be [23 its, 31 its, 1222 its, Did not converge] in that order.
summary %>% 
  kbl(caption = "Summary on the corrupted mat") %>% kable_classic(full_width = F, html_font = "Cambria", position = "center") %>% 
  kable_styling(bootstrap_options = c("hover", "condensed"), fixed_thead = T)
