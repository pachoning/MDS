library(profvis)
source("final_methods/divide_conquer_mds.R")
source("final_methods/fast_mds.R")
source("final_methods/gower_interpolation_mds.R")

iterations = 3
s = 10
n = s^(iterations + 2)
cols = 4
x = matrix(data = rnorm(n*cols), ncol = cols)
l = (n*s^iterations)^(1/(iterations + 1))
k = 4
tie = s
dist_fn = stats::dist

profvis({
  divide_conquer_mds(x, l, tie, k)
})


profvis({
  fast_mds(x, l, s, k)
})


profvis({
  gower_interpolation_mds(x, l, k)
})
