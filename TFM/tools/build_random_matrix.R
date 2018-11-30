build_matrix <- function(
  n,
  p,
  corr_coef
){
  S <- matrix(corr_coef , ncol = p, nrow = p)
  diag(S) = 1
  r <- eigen(S)
  V <- r$vectors
  lam <- r$values
  A <- V %*% diag( lam^(1/2) ) %*% t(V)  
  
  return( matrix( rnorm( n*p ), nrow = n ) %*% A )
  
}

