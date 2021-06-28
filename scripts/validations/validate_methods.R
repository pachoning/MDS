source("tools/load_libraries.R")
source("tools/procrustes.R")
source("final_methods/classical_mds.R")

source("final_methods/divide_conquer_mds.R")
source("final_methods/divide_conquer_small_mds.R")
source("final_methods/fast_mds.R")

source("final_methods/gower_interpolation_mds.R")
source("final_methods/gower_interpolation_rowwise_mds.R")
source("final_methods/gower_interpolation_kgropus_mds.R")
source("final_methods/gower_interpolation_dist_mds.R")
source("final_methods/gower_interpolation_big_mds.R")

# Generate data
n_rows <- 1000
var_vector <- c(5, 5, 1)
mean_vector <- c(0, 0, 0)
n_cols <- length(var_vector)
x <- mapply(rnorm, n = n_rows, mean_vector) %*% diag(sqrt(var_vector))
dim(x)
var(x)
cor(x)

s <- 2*n_cols
tie <- s
k <- n_cols
l <- 200
dist_fn = stats::dist

########################################################################################################################
##### Results 
algorithm <- gower_interpolation_big_mds
results <- algorithm(x = x,l = l, tie = tie, k = k, dist_fn = dist_fn, s = s)
procrustes <- perform_procrustes(x = results$points,
                                 target = x, 
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
