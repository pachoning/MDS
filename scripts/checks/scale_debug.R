source("final_methods/classical_mds.R")
source("final_methods/interpolation_mds.R")
source("final_methods/landmark_mds.R")
source("final_methods/fast_mds.R")
source("final_methods/pivot_mds.R")
source("final_methods/reduced_mds.R")
source("final_methods/divide_conquer_mds.R")

n_col <- 5
n_row <- 1000
x <- matrix(
  data = rnorm(n_col * n_row),
  nrow = n_row
)

l_val <- 100
x <- x %*% diag(c(15, 10, 1, 1, 1)^.5)
#x <- x %*% (matrix(1,nrow = 5, ncol = 5) + diag(c(15,10,1,1,1)^.5))
var(x)


results_classical <- classical_mds(x, k = 2)
resuls_interpolation <- interpolation_mds(x, l = l_val, k = 2, n_cores = 1)
results_landmark <- landmark_mds(x, ndim = 2, num_landmarks = l_val)
results_fast <- fast_mds(x, l = l_val, k = 2, n_cores = 1, s = 10)
results_pivot <- pivot_mds(x, pivots = l_val, ndim = 2)
results_reduced <- reduced_mds(x, l = l_val, k = 2)
results_dc <- divide_conquer_mds(x, l =l_val, tie = 10, k = 2, n_cores = 1)

#plot(results_classical$points, xlim = c(-15, 15), ylim = c(-10, 10))
idex_sample <- sample(n_row, n_row)
op <- par(mfrow = c(2, 3))
plot(resuls_interpolation$points[idex_sample, ], xlim = c(-15, 15), ylim = c(-10, 10))
plot(results_landmark$points[idex_sample, ], xlim = c(-15, 15), ylim = c(-10, 10))
plot(results_fast$points[idex_sample, ], xlim = c(-15, 15), ylim = c(-10, 10))
plot(results_pivot$points[idex_sample, ], xlim = c(-15, 15), ylim = c(-10, 10))
plot(results_reduced$points[idex_sample, ], xlim = c(-15, 15), ylim = c(-10, 10))
plot(results_dc$points[idex_sample, ], xlim = c(-15, 15), ylim = c(-10, 10))
par(op)
plot(results_pivot$points)
