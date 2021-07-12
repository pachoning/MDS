source("tools/load_libraries.R")
source("tools/procrustes.R")
source("final_methods/classical_mds.R")
source("final_methods/divide_conquer_mds.R")
source("final_methods/gower_interpolation_mds.R")
source("final_methods/fast_mds.R")

# Generate data
n_rows <- 10^5
n_cols <- 100
n_dominant_dimension <- 5
var_vector <- c(rep(15, n_dominant_dimension), rep(1, n_cols - n_dominant_dimension))
mean_vector <- rep(0, n_cols)
x <- mapply(rnorm, n = n_rows, mean_vector) %*% diag(sqrt(var_vector))
dim(x)
#var(x)
#cor(x)

s <- 2*n_cols
tie <- s
k <- n_cols
l <- 1500
n_row_partition <- 1
dist_fn <-  stats::dist

##############################################################################################################
##### Results 
algorithm <- gower_interpolation_mds
all_n_row_partition <- c(1, seq(from = 100, to = 1500, by = 100), (n_rows - l))

for (n_row_partition in all_n_row_partition) {
  ini_time <- proc.time()
  results <- algorithm(x = x, l = l, tie = tie, k = k, dist_fn = dist_fn, s = s, n_row_partition = n_row_partition)
  final_time <- proc.time()
  elapsed_time <- final_time - ini_time
  elapsed_time <- elapsed_time[3]
  msg <- paste0("For value: ", n_row_partition, " the time has been: ", elapsed_time)
  message(msg)
}

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
