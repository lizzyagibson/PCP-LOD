\name{angeldust}
\alias{angeldust}
\title{Principal Component Pursuit}
\usage{
pcp(D, mu)
}
\description{
PCP decomposes an exposure matrix into a low rank matrix of exposure patterns and a sparse matrix of unique or extreme exposure events.

L = Low rank matriz,
S = Sparse matrix,
Lambda = lambda, penalization parameter set to 1/sqrt(nrows of original data matrix),
Mu = mu, penalization parameter chosen by user,
obj_value = objective function value,
error = D - L - S, lost information
}

\examples{
pcp(D, mu)
}
