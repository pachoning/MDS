source("final_methods/fast_mds.R")

x <- matrix(rnorm(10*8), ncol = 8)

l <- 7
s <- 1
k <- 1
fast_mds(x, l, s = s, k = k, n_cores = 1)
