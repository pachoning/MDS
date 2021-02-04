source("tools/classical_mds.R")
source("tools/fast_mds.R")
source("tools/divide_conquer_mds.R")
source("tools/gower_interpolation_mds.R")
source("tools/procrustes.R")

# Generate data
n_rows <- 10000
var_vector <- c(5, 5, 1)
n_cols <- length(var_vector)
x <- matrix(rnorm(n_rows*n_cols), nrow = n_rows) %*% diag(var_vector)
dim(x)
var(x)
cor(x)

s <- 2*n_cols
k <- n_cols
l <- 100

##### Results for divide and conquer mds
divide_results <- divide_conquer_mds(x = x,l = l, tie = s, k = k, dist_fn = stats::dist)
divide_proc <- perform_procrustes(x = divide_results$points, target = x, matrix_to_transform = divide_results$points, 
                                 translation = FALSE, dilation = FALSE)

# Correlation of principal coordinates
cor(divide_proc[,1], x[, 1])

# Eigenvalues
sd(divide_results$points[, 1])
sqrt(divide_results$eigen)



##### Results for fast mds
fast_results <- fast_mds(x = x,l = l, s = s, k = k, dist_fn = stats::dist)
fast_proc <- perform_procrustes(x =fast_results$points, target = x, matrix_to_transform = fast_results$points, 
                                 translation = FALSE, dilation = FALSE)

# Correlation of principal coordinates
cor(fast_proc[,1], x[, 1])

# Eigenvalues
sd(fast_results$points[, 1])
sqrt(fast_results$eigen)

##### Results for Gower mds
gower_results <- gower_interpolation_mds(x = x, l = l, k = k, dist_fn = stats::dist)
gower_proc <- perform_procrustes(x = gower_results$points, target = x, matrix_to_transform = gower_results$points, 
                                translation = FALSE, dilation = FALSE)

# Correlation of principal coordinates
cor(gower_proc[,1], x[, 1])

# Eigenvalues
sd(gower_results$points[, 1])
sqrt(gower_results$eigen)

