library(profvis)
source("final_methods/divide_conquer_mds.R")
source("final_methods/fast_mds.R")
source("final_methods/gower_interpolation_mds.R")

iterations <- 3
k <- 5
s <- 2*k
#n <- s^(iterations + 2)
n <- 10^5
cols <- 100
x <- matrix(data = rnorm(n*cols), ncol = cols)
#l <- (n*s^iterations)^(1/(iterations + 1))
l <- 1500
tie <- s
dist_fn <- stats::dist

profvis({
  divide_conquer_mds(x, l, tie, k)
})


profvis({
  fast_mds(x, l, s, k)
})


profvis({
  gower_interpolation_mds(x, l, k)
})

profvis({
  gower_interpolation_rowwise_mds(x, l, k)
})

