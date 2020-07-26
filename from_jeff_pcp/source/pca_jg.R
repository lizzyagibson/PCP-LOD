## do pca with given rank and return a fitted matrix
## i'm making some choices here, especially when the matrix has missing values
## i don't know if those are good choices though

pca_jg = function(mat, rank = 1) {
  
  I = dim(mat)[1]
  D = dim(mat)[2]
  
  data = replace(mat, mat == -1, NA)
  
  mean_vec = apply(data[complete.cases(data),], 2, mean, na.rm = TRUE)
  data.tilde = data - matrix(mean_vec, I, D, byrow = TRUE)
  
  eigen_decomp =
    cov(data.tilde, use = "complete.obs") %>% 
    eigen
  
  est_patterns = eigen_decomp$vectors[, 1:rank]
  
  est_scores = matrix(NA, I, rank)
  est_mat = matrix(NA, I, D)
  
  for(i in 1:I){
    est_scores[i, ] = lm(data.tilde[i, ] ~ 0 + est_patterns) %>% coef
    est_mat[i, ] = mean_vec + est_scores[i, ] %*% t(est_patterns)
  }
  
  est_mat
    
}

#############

pca_lz <- function(mat, rank = 1) {
  
  mean_vec = apply(data, 2, mean, na.rm = TRUE)
  
  pca <- prcomp(data)
  
  pca_pred  <- as.matrix(pca$x[,1:rank]) %*% t(pca$rotation)[1:rank,] + matrix(rep(mean_vec,each=100),nrow=100)
  
  pca_pred
}
