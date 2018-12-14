x = matrix(rnorm(100*5), ncol = 5)
nrow_x = nrow(x)
ncol_x = ncol(x)


distance_matrix = 
  In = diag(nrow_x)
ones_vector = rep(1, nrow_x)
J = In - 1/nrow_x*ones_vector %*% t(ones_vector)
delta_matrix = 