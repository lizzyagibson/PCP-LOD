# Author: Lawrence Chillrud (lgc2139@cumc.columbia.edu)
# Date: 3/10/21

# Contents:
# 0. PACKAGE IMPORTS
# 1. SIMULATE DATA MATRICES
# 2. SET UP GRIDSEARCHES
# 3. CONDUCTING GRIDSEARCHES
  # a. A VANILLA GRID SEARCH
  # b. A SEARCH OF JUST LAMBDAS
  # c. A BAYES SEARCH WITH NON-CONVEX PCP
  # d. A RANDOM SEARCH WITH NON-CONVEX PCP
# 4. CHOOSING OPTIMAL PARAMS AND RUNNING NON-CONVEX PCP

######## 0. PACKAGE IMPORTS ########
library(PCPhelpers)
library(pcpr)
library(magrittr)

######## 1. SIMULATE DATA MATRICES ########

# Here, we simulate a 50x10 rank3 matrix (with no sparse component, some noise (stdev = .1))
# named "mat" (line 13)
n <- 50
p <- 10
data <- sim_data(1, n, p, 3, sigma=0.1) # PCPhelper function that returns a list = M, L, S, and Z
mat <- data$M

# Here, we subject mat to an artificial LOD yielding "mat.min1"
# mat.min1 has entries < LOD encoded as -1
# The LOD for each column is recorded in "delta"
delta <- apply(mat, 2, quantile, probs = .25)
mat.min1 <- corrupt_mat(mat, cols = 1:ncol(mat), limit = .25, fill = "-1")

######## 2. SET UP GRIDSEARCHES ########
lambda <- seq(.1, .5, length.out = 5)
mu <- seq(.5, 4, length.out = 5)
rank <- 1:5

grid.l <- data.frame(lambda = lambda) # a grid only searching through lambdas
grid.lm <- expand.grid(lambda = lambda, mu = mu) # a grid searching through both lambda and mu
grid.lr <- expand.grid(lambda = lambda, r = rank) # a grid searching through only lambda and rank
grid.lmr <- expand.grid(lambda = lambda, mu = mu, r = rank) # a grid searching through lambda, mu and rank (for non-convex PCP)

######## 3. CONDUCTING GRIDSEARCHES ########

#### 3a. A VANILLA GRID SEARCH ####

# example with a vanilla grid search:
# here I'm using the basic grid_search function that naively / brute force looks at the entire grid_df
# we're using root_pcp_na_nonnegL, aka Root PCP that can handle NAs and applies a non-negativity constraint on L
# we're searching through different combinations of lambda and mu with the grid.lm data frame we made.
vanilla_search <- grid_search_cv(mat = mat, pcp_func = root_pcp_na_nonnegL, grid_df = grid.lm,
                                 cores = 4, runs = 5)

print_gs(vanilla_search$formatted) # viewing the main results
View(vanilla_search$raw) # notice that the raw data frame also has information on the rank of L and sparsity of S for each run
vanilla_search$constants # this will be an empty list, since we didn't supply any constant parameters for PCP.

#### 3b. A SEARCH OF JUST LAMBDAS ####

# this example is identical to the one above, except: 
# I'm only looking through values of lambda (using grid.l), while holding mu at a constant value of 3.
lambdas_search <- grid_search_cv(mat = mat, pcp_func = root_pcp_na_nonnegL, grid_df = grid.l,
                                 cores = 4, runs = 5, mu = 3) # <- mu is held at a constant 3 here.

View(lambdas_search$formatted) # viewing the main results (we can no longer use print_gs, since at the moment it's just for {lambda x mu} grid searches)
View(lambdas_search$raw) # the raw data frame also has information on the rank of L and sparsity of S for each run
lambdas_search$constants # Now this will be a list containing mu = 3, reminding us that we held mu to a constant 3 throughout the search.

#### 3c. A BAYES SEARCH WITH NON-CONVEX PCP ####

# in this example, we'll use non-convex PCP, which takes the desired rank of the recovered L matrix as a parameter, via the "r" argument
# we will conduct a bayes "intelligent" search of our parameter space {lambda x mu x r}
# therefore, we'll pass the grid.lmr data frame we made above. 
# we'll also use the LOD version of non-convex PCP (along with NA values / non-negative L),
# so we need to supply the LOD argument that non-convex PCP will use, but hold it constant, like we did mu 
# in the example above. We also need to supply mat.min1, instead of mat, which has LOD information encoded
# as -1 values.
sample_search <- bayes_search_cv(mat = mat.min1, pcp_func = root_pcp_noncvx_nonnegL_na_lod, grid_df = grid.lmr,
                                 init_evals = 10, bayes_evals = 6, 
                                 cores = 4, runs = 5, LOD = delta) # <- here we hold LOD constant

# we can even continue the search by passing bayes the formatted results from the last one:
noncvx_search2 <- bayes_search_cv(mat = mat.min1, pcp_func = root_pcp_noncvx_nonnegL_na_lod, grid_df = sample_search$formatted,
                                 init_evals = 10, bayes_evals = 1, 
                                 cores = 4, runs = 5, LOD = delta)
print_gs(noncvx_search2$formatted) # viewing the main results
class(noncvx_search2$formatted)
#### 3d. A RANDOM SEARCH WITH NON-CONVEX PCP ####

# We can hold multiple PCP parameters constant by just passing them to our search_cv 
# function instead of giving them their own column in grid_df. 
# For example, let's hold LOD and mu constant, and only look through lambda and the rank:
random_search <- random_search_cv(mat = mat.min1, pcp_func = root_pcp_noncvx_nonnegL_na_lod, grid_df = grid.lr,
                                 n_evals = 7, cores = 4, runs = 5, 
                                 mu = 3, LOD = delta) # <- here we hold both mu and LOD constant

######## 4. CHOOSING OPTIMAL PARAMS AND RUNNING NON-CONVEX PCP ########

# in practice, we choose our "optimal" parameters by examining the L.rank and S.sparsity metrics from the
# raw dataframe outputted by the search, finding runs where those metrics make intuitive sense for the application.
# then we narrow down parameters by looking at those parameter settings where "value" in the formatted dataframe is lowest.

# for these toy examples that don't have any meaning, we'll just find our optimal parameters by taking those that have
# the minimum "value" (aka relative recovery error) from the formatted dataframe.

best.idx <- which.min(noncvx_search$formatted$value)
optimal.lambda <- noncvx_search$formatted$lambda[best.idx]
optimal.mu <- noncvx_search$formatted$mu[best.idx]
optimal.r <- noncvx_search$formatted$r[best.idx]

noncvx.out <- root_pcp_noncvx_nonnegL_na_lod(mat.min1, lambda = optimal.lambda, mu = optimal.mu, r = optimal.r, 
                                             LOD = delta, verbose = T)

View(noncvx.out$L)
View(noncvx.out$S)