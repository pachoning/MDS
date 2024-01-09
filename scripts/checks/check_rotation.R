source("final_methods/landmark_mds.R")
source("final_methods/interpolation_mds.R")
source("final_methods/fast_mds.R")
source("final_methods/divide_conquer_mds.R")
source("final_methods/pivot_mds.R")
source("final_methods/reduced_mds.R")


n_rows <- 1000
n_col <- 10
k <- 3
l <- 100
X <- matrix(rnorm(n_rows * n_col), ncol = n_col)
debug(fast_mds)
results <- fast_mds(X, l=500, s = 20, k = k, n_cores = 1)
cor(landmark_results$points)

debug(reduced_mds)
results <- reduced_mds(X, k = k, l =l)
