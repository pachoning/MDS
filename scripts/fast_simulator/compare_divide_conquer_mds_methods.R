source("tools/divide_conquer_mds.R")
source("tools/divide_conquer_mds_candidate_1.R")
source("tools/fast_mds.R")

n_sim <- 20
num_rows <- 100000
num_cols <- 100
k <- 4
l <- 200
tie <- 10
dist_fn <- stats::dist
time_divide_conquer_mds <- c()
time_divide_conquer_mds_efficient <- c()

for (i_sim in 1:n_sim) {
  x <- matrix(rnorm(num_rows*num_cols), ncol = num_cols)
  message(paste0("Working on simulation: ", i_sim))
  init_time <- proc.time()
  mds_efficient <- divide_conquer_mds_efficient(x = x, l = l, tie = tie, k = k, dist_fn = stats::dist)
  elapsed_time <- proc.time() - init_time
  time_divide_conquer_mds_efficient <- c(time_divide_conquer_mds_efficient, elapsed_time[3])

  init_time <- proc.time()
  mds_efficient <- divide_conquer_mds(x = x, l = l, tie = tie, k = k, dist_fn = stats::dist)
  elapsed_time <- proc.time() - init_time
  time_divide_conquer_mds <- c(time_divide_conquer_mds, elapsed_time[3])
}

mean(time_divide_conquer_mds)
mean(time_divide_conquer_mds_efficient)




x <- matrix(rnorm(num_rows*num_cols), ncol = num_cols)
init_time <- proc.time()
res <- divide_conquer_mds(x = x, l = l, tie = tie, k = k, dist_fn = stats::dist)
elapsed_time <- proc.time() - init_time
message(paste0("------------ Total time: ", elapsed_time[3]))

init_time <- proc.time()
res <- fast_mds(x = x, l = l, s = tie, k, dist_fn = stats::dist)
elapsed_time <- proc.time() - init_time
message(paste0("------------ Total time: ", elapsed_time[3]))


