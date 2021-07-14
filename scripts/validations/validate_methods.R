source("tools/load_libraries.R")
source("tools/procrustes.R")
source("final_methods/classical_mds.R")
source("final_methods/divide_conquer_mds.R")
source("final_methods/gower_interpolation_mds.R")
source("final_methods/fast_mds.R")

# Generate data
n_rows <- 10^6
n_cols <- 100
n_dominant_dimension <- 10
var_vector <- c(rep(15, n_dominant_dimension), rep(1, n_cols - n_dominant_dimension))
mean_vector <- rep(0, n_cols)
#x <- mapply(rnorm, n = n_rows, mean_vector) %*% diag(sqrt(var_vector))
x <- matrix(data = rnorm(n_rows * n_cols), nrow = n_rows) %*% diag(sqrt(var_vector))
dim(x)
#var(x)
#cor(x)

s <- 2*n_dominant_dimension
tie <- s
k <- n_dominant_dimension
l <- 200
dist_fn <-  stats::dist
n_cores <- 5

##############################################################################################################
##### Results 
algorithm <- gower_interpolation_mds
results <- algorithm(x = x, l = l, tie = tie, k = k, dist_fn = dist_fn, s = s, n_cores = n_cores)
procrustes <- perform_procrustes(x = results$points,
                                 target = x[, 1:n_dominant_dimension, drop = FALSE], 
                                 matrix_to_transform = results$points,
                                 translation = FALSE)

# Checking it is a MDS solution
nrow(results$points)
var(results$points)
cor(results$points)
apply(results$points, 2, mean)

# Correlation of principal coordinates
cor(procrustes[,1], x[, 1])

# Eigenvalues
var(results$points)
results$eigen

# Cheching distance is preserved
dist(x[1:5, ])
dist(results$points[1:5, ])

# Checking GOF
results$GOF
