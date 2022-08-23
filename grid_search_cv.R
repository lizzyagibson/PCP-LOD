# source pcp function
source("pcpr.R")

# Code for grid search, with necessary nested functions
# Rough way to determine the sparsity of a given matrix.
sparsity <- function(mat, tol = .00001) {
  mat[abs(mat) < tol] <- 0
  100 - (100 * sum((mat != 0) * 1) / prod(dim(mat)))
}

# Corrupts a data matrix
corrupt_mat <- function(mat, cols, limit, fill="NA") {
  mat[, cols] <- apply(mat[, cols, drop=FALSE], FUN=impute, limit=limit, fill=fill, MARGIN=2)
  mat
}

# Corrupts a data matrix by knocking out a random specified percentage of entries.
corrupt_mat_randomly <- function(mat, k, perc_b) {
  nvals_to_corrupt <- floor(perc_b*prod(dim(mat)))
  mat.vec <- as.vector(mat)
  mask <- rep(0, length(mat.vec))
  
  pool <- which(mat.vec >= 0)
  if (nvals_to_corrupt > length(pool)) {
    corrupted <- pool
  } else {
    set.seed(k)
    corrupted <- sample(pool, nvals_to_corrupt, replace = FALSE)
  }
  
  mask[corrupted] <- 1
  mat.vec[corrupted] <- NA
  
  rows <- nrow(mat)
  cols <- ncol(mat)
  
  ret.mat <- matrix(mat.vec, nrow = rows, ncol = cols)
  ret.mask <- matrix(mask, nrow = rows, ncol = cols)
  
  list(cor.mat = ret.mat, cor.mask = ret.mask)
}

# Evaluates a given setting of lambda and mu on a given matrix with a given PCP function. 
eval_params <- function(
    seed,
    mat,
    pcp_func,
    perc_b,
    eval_params_index,
    ...
) {
  
  pcp_args <- list(...)
  
  cor_mat <- corrupt_mat_randomly(mat, k = seed, perc_b = perc_b) # a list containing cor.mat and cor.mask
  mask <- cor_mat$cor.mask
  pcp_out <- do.call(pcp_func, c(pcp_args, list(D = cor_mat$cor.mat, verbose = F)))
  mat[is.na(mat)] <- 0
  score <- norm((mat - pcp_out$L - pcp_out$S)*mask, "F") / norm(mat * mask, "F")
  L.rank <- Matrix::rankMatrix(pcp_out$L, tol = 1e-04)
  S.sparsity <- sparsity(pcp_out$S, tol = 1e-04)
  its <- pcp_out$final_iter
  
  c(as.numeric(pcp_args[eval_params_index]), seed, score, L.rank, S.sparsity, its)
}

# Conducts a cross-validated grid search of the parameters for Principle Component Pursuit (PCP).
grid_search_cv <- function(
    mat, 
    pcp_func,
    grid_df, 
    cores = NULL,
    perc_b = 0.2,
    runs = 100,
    progress_bar = TRUE,
    file = NULL,
    ...) 
{
  
  # setting up the parallel programming
  if (is.null(cores)) {
    cores <- parallel::detectCores()
  } 
  cl <- snow::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  
  # initialization
  metrics <- c("value", "L.rank", "S.sparsity", "iterations")
  for (metric in metrics) {
    if (!metric %in% names(grid_df)) {
      grid_df[, metric] <- as.numeric(NA)
    }
  }
  
  param_names <- grid_df %>% dplyr::select(!tidyselect::all_of(metrics)) %>% colnames()
  points_to_eval <- which(is.na(grid_df$value))
  params <- data.frame(grid_df[points_to_eval, param_names], rep(1:runs, each = length(points_to_eval)), row.names = NULL)
  colnames(params) <- c(param_names, "seed")
  constant_params <- list(...)
  
  # progress bar setup
  if (progress_bar) {
    pb <- txtProgressBar(min=0, max=nrow(params), width=50, style=3)
    progress <- function(p) setTxtProgressBar(pb, p)
    opts <- list(progress=progress)
  } else {
    opts <- list()
  }
  
  # grid search
  cv <- foreach(i = iterators::icount(nrow(params)), .options.snow=opts, .combine = cbind,
                .inorder = FALSE, 
                .export=c("eval_params", "corrupt_mat_randomly", "prox_nuclear", "prox_l1", "proj_rank_r",
                          "prox_fro", "sparsity")) %dopar% {
    do.call(what = 'eval_params', c(as.list(params[i,]), constant_params, list(mat = mat, pcp_func = pcp_func, perc_b = perc_b, eval_params_index = param_names)))
  }
  
  # close the progress bar and stop cluster
  if (progress_bar) close(pb)
  snow::stopCluster(cl)
  
  # format the output
  rownames(cv) <- c(colnames(params), metrics)
  cv <- as.data.frame(t(cv))
  
  cv.formatted <- cv %>% 
    dplyr::group_by_at(all_of(param_names)) %>% 
    dplyr::summarise(value = mean(value), L.rank = mean(L.rank), S.sparsity = mean(S.sparsity), iterations = mean(iterations)) %>% 
    tidyr::unite(param_setting, all_of(param_names))
  
  grid_df.formatted <- grid_df %>% 
    tidyr::unite(param_setting, all_of(param_names))
  
  grid_df.formatted[match(cv.formatted$param_setting, grid_df.formatted$param_setting), ] <- cv.formatted
  
  grid_df <- grid_df.formatted %>% 
    tidyr::separate(param_setting, into = param_names, sep = "_", convert = TRUE)
  
  # save
  if (!is.null(file)) save(cv, grid_df, file = file)
  
  # return
  list(raw = cv, formatted = grid_df, constants = constant_params)
}
