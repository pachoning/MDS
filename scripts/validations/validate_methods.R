source("tools/load_libraries.R")
source("tools/procrustes.R")
source("final_methods/classical_mds.R")
source("final_methods/fast_mds.R")
source("final_methods/divide_conquer_mds.R")
source("final_methods/gower_interpolation_mds.R")


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

########################################################################################################################
##### Results for divide and conquer mds
divide_results <- divide_conquer_mds(x = x,l = l, tie = tie, k = k, dist_fn = stats::dist)
divide_proc <- perform_procrustes(x = divide_results$points, target = x, matrix_to_transform = divide_results$points, 
                                 translation = FALSE)

# Checking it is a MDS solution
var(divide_results$points)
cor(divide_results$points)
apply(divide_results$points, 2, mean)

# Correlation of principal coordinates
cor(divide_proc[,1], x[, 1])

# Eigenvalues
var(divide_results$points)
divide_results$eigen

# Cheching distance is preserved
dist(x[1:5, ])
dist(divide_results$points[1:5, ])

# Checking GOF
divide_results$GOF
########################################################################################################################
##### Results for fast mds
fast_results <- fast_mds(x = x,l = l, s = s, k = k, dist_fn = stats::dist)
fast_proc <- perform_procrustes(x = fast_results$points, target = x, matrix_to_transform = fast_results$points, 
                                 translation = FALSE)

# Checking it is a MDS solution
cov(fast_results$points)
cor(fast_results$points)
apply(fast_results$points, 2, mean)

# Correlation of principal coordinates
cor(fast_proc[,1], x[, 1])

# Eigenvalues
var(fast_results$points)
fast_results$eigen

# Cheching distance is preserved
dist(x[1:5, ])
dist(fast_results$points[1:5, ])

# Checking GOF
fast_results$GOF
########################################################################################################################
##### Results for Gower mds
gower_results <- gower_interpolation_mds(x = x, l = l, k = k, dist_fn = stats::dist)
gower_proc <- perform_procrustes(x = gower_results$points, target = x, matrix_to_transform = gower_results$points, 
                                translation = FALSE)

# Checking it is a MDS solution
cov(gower_results$points)
cor(gower_results$points)
apply(gower_results$points, 2, mean)


# Correlation of principal coordinates
cor(gower_proc[,1], x[, 1])

# Eigenvalues
var(gower_results$points)
gower_results$eigen

# Cheching distance is preserved
dist(x[1:5, ])
dist(gower_results$points[1:5, ])

# Checking GOF
gower_results$GOF
