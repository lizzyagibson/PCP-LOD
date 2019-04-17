# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

########################################################################
########################################################################

soft_thresholding <- function(v, lambda) {
  myzero <- matrix(data = 0, ncol = ncol(v), nrow = nrow(v))
  w <- sign(v) * pmax(abs(v) - lambda, myzero)
  # If absolute value is less than set lambda value, push to zero
  # If absolute value is greater than lambda, new value is difference
  # W is either zero or the difference between v and lambda
  w
}

soft_thresholding_diag <- function(v, lambda) {
  myzero <- vector("numeric", length = length(v))
  w <- sign(v) * pmax(abs(v) - lambda, myzero)
  # If absolute value is less than set lambda value, push to zero
  # If absolute value is greater than lambda, new value is difference
  # W is either zero or the difference between v and lambda
  w
}

########################################################################
########################################################################

singular_value_threshold <- function(M, lambda) {

  USV <- svd(M)
  # Break SVD into separate matrices
  U <- USV$u
  # U is each persons experience of the source
  # Matrix whose columns contain the left singular vectors of M. Dimension c(n, nu).
  sv <- USV$d
  # Diagonal matrix of singular values, sorted decreasingly.
  V <- USV$v
  # Exposures and sources contribution

  N <- U %*% diag(soft_thresholding_diag(sv, lambda)) %*% t(V)
  # Create a new version of the input matrix by putting the SVD back together
  # With new singular value diagonal matrix
  # Singular values greater than assigned lambda are pushed to zero
  # Singular value is relative scaling of how much a source is contributing
  # Zero singular values will make this new matrix lower rank

  v  <- sum(soft_thresholding_diag(sv, lambda))
  # Sum the singular values, with those less than assigned lambda pushed to zero

  svt <- list(N = N, v = v)
  # Output new data matrix and sum of singular values

  svt
}

########################################################################
########################################################################

pcp <- function(D, mu) {

  D <- as.matrix(D)
  # Dataframe needs to be a matrix

  m <- nrow(D)
  n <- ncol(D)
  lambda <- 1/sqrt(m)

  S <- matrix(0, nrow = m, ncol = n)
  L <- matrix(0, nrow = m, ncol = n)
  # First iteration starts with empty matrices same size as original data matrix

  iter <- 0
  done <- FALSE
  MAX_ITER <- 100
  # Maximum number of iterations

  while (!done) {

    iter <- iter + 1
    #Loop through this function, update every time

    svt <- singular_value_threshold((D - S), 1/mu)
    # First iteration is on original data matrix
    # Following iterations are on D - S = Low rank matrix
    # Singular values less than 1/mu are pushed to zero
    # Singular values are either zero or the remainder after subtracting 1/mu

    L <- svt[[1]] #svt$N
    v <- svt[[2]]
    # Outputs thresholded Low rank matrix and sum of thresholded singular values (v)

    S <- soft_thresholding((D - L), lambda/mu)
    # First iteration is on original data matrix
    # Following iterations are on D - L = Sparse matrix
    # All values in Sparse matrix less than lambda/mu are pushed to zero
    # Values are either zero or the remainder after subtracting lambda/mu

    error <- D - L - S
    # D - L - S should be close to zero, error
    # Parts of exposure matrix we are losing

    obj <- v + lambda * sum(abs(S)) + (mu/2) * norm((D - L - S), type = "F")^2
    # Objective function will decrease with each iteration
    # Larger error matrix => larger objective function
      # Penalizing error in prediction
    # v is sum of singular values from Low rank matrix
      # Penalizing larger values
    # abs(S) means less sparse sparse matrix is more penalized

    print(paste(iter, "Obj:", obj))

    if (iter >= MAX_ITER) {done <- TRUE}

  }
  list(L = L, S = S, Lambda = lambda, Mu = mu, obj_value = obj, error = error)
}

