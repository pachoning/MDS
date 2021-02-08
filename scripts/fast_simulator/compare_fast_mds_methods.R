source("tools/fast_mds.R")
source("tools/fast_mds_candidate_1.R")

n_sim <- 20
num_rows <- 100000
num_cols <- 10
k <- 4
l <- 200
s <- 10
dist_fn <- stats::dist
time_fast <- c()
time_candidate_1 <- c()

for (i_sim in 1:n_sim) {
  x <- matrix(rnorm(num_rows*num_cols), ncol = num_cols)
  message(paste0("Working on simulation: ", i_sim))
  init_time <- proc.time()
  fast <- fast_mds(x = x, l = l, s = s, k = k, dist_fn = stats::dist)
  elapsed_time <- proc.time() - init_time
  time_fast <- c(time_fast, elapsed_time[3])
  
  init_time <- proc.time()
  mds_efficient <- fast_mds_candidate_1(x = x, l = l, s = s, k = k, dist_fn = stats::dist)
  elapsed_time <- proc.time() - init_time
  time_candidate_1 <- c(time_candidate_1, elapsed_time[3])
}

mean(time_fast)
mean(time_candidate_1)
