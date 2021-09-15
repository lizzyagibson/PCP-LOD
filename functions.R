# packages used
packages <- c( "haven", "tidyverse", "RColorBrewer", "LearnBayes", "janitor", "broom", 
               "ggrepel","ggpubr", "FactoMineR", "osqp", "CVXR", "factoextra",
               "GGally", "textshape",  "grid", "pracma", 
               "heatmaply","ggsci", "patchwork", "here",
               "snow", "doSNOW", "GPfit", "kableExtra")

# if these aren't installed, install them
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {install.packages(new.packages)}

# load all packages
lapply(packages, library, character.only = TRUE)

# Change this absolute path to your own
# install.packages("/Users/lizzy/pcpr", repos = NULL, type="source")
library(pcpr)

# Change this absolute path to your own
# install.packages("/Users/lizzy/experiments/law_experiments/PCPhelpers", repos = NULL, type="source")
library(PCPhelpers)

theme_set(theme_bw(base_size = 20) + theme(legend.position = "bottom",
                                           strip.background =element_rect(fill="white")))

# L2 relative error
# WITH error handling
# if the matrices are not the same size, NA
get_relerror <- function(x,y) {
  
  if(any(is.na(x)) | any(is.na(y))) {return(NA)} else{
    
    x = as.matrix(x)
    y = as.matrix(y)
    
    # patterns should be columns
    if (nrow(x) < ncol(x)) {x <- t(x)}
    if (nrow(y) < ncol(y)) {y <- t(y)}  
    
    if(ncol(x) != ncol(y)) { return(NA) } else{
      return(norm(x-y, "F")/norm(x, "F"))
    }
  }
}

# Factor correspondence
# WITH error handling
# if matrices are not same size, NA
# If nn is false, we allow *signed permutations*
factor_correspondence <- function (A, B, nn = FALSE) {
  # This is all from the CVXR package, which is kind of its own language
  # for convex optimization
  G <- t(B) %*% A
  n <- nrow(G)
  
  # Step 1. Define the variable to be estimated
  # Pi -- the permutation or signed permutation matrix
  Pi <- Variable(n,n)
  
  # Step 2. Define the objective to be optimized
  objective <- Maximize(base::sum(Pi * G))
  
  if (nn) {
    # Step 2.5. Subject to these constraints
    constX = list()
    for (i in 1:nrow(G)) {
      constX <-  c(constX, list(base::sum(Pi[,i]) == 1))
      constX <-  c(constX, list(base::sum(Pi[i,]) == 1))
    }
    constX <- c(constX, Pi >= 0)
  } else {
    # % allow sign flips 
    # Step 3. vector l1 norms along rows and columns
    constX = list()
    for (i in 1:nrow(G)) {
      constX <-  c(constX, list(base::sum(abs(Pi[,i])) <= 1))
      constX <-  c(constX, list(base::sum(abs(Pi[i,])) <= 1))
    }
  }  
  # Step 3. Create a problem to solve
  problem <- Problem(objective, constraints = constX)
  
  # Step 4. Solve it!
  result <- tryCatch({
    solve(problem)
  }, warning = function(warning_condition) {
    message('Caught a warning!')
  }, error = function(error_condition) {
    message('Caught an error!')
  })
  
  # Step 5. Extract solution and objective value
  perm = tryCatch({
    round(result$getValue(Pi), 0)
  }, warning = function(warning_condition) {
    message('Caught a warning!')
    diag(ncol(B))
  }, error = function(error_condition) {
    message('Caught an error!')
    diag(ncol(B))
  }, finally={
    diag(ncol(B))
  })
  
  e <- norm(B,'f')^2 + norm(A,'f')^2 - 2 * base::sum(perm * G)
  # e -- the sum of squared errors under the best calibration \Pi
  
  # New matrix with best order
  newB <- B %*% perm
  
  return(newB)
}

# Get triangle for correlation heatmap
get_lower_tri <-function(x){
  x[upper.tri(x)] <- NA
  return(x)
}
