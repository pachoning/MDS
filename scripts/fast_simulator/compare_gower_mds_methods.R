source("tools/gower_interpolation_mds.R")
source("tools/gower_interpolation_mds_candidate_1.R")

n_sim <- 20
num_rows <- 100000
num_cols <- 10
k <- 4
l <- 200
dist_fn <- stats::dist
time_gower_mds <- c()
time_gower_candidate_1 <- c()

for (i_sim in 1:n_sim) {
  x <- matrix(rnorm(num_rows*num_cols), ncol = num_cols)
  message(paste0("Working on simulation: ", i_sim))
  init_time <- proc.time()
  gow <- gower_interpolation_mds(x = x, l = l, k = k, dist_fn = stats::dist)
  elapsed_time <- proc.time() - init_time
  time_gower_mds <- c(time_gower_mds, elapsed_time[3])
  
  init_time <- proc.time()
  gower_cand_1 <- gower_interpolation_mds_candidate_1(x = x, l = l, k = k, dist_fn = stats::dist)
  elapsed_time <- proc.time() - init_time
  time_gower_candidate_1 <- c(time_gower_candidate_1, elapsed_time[3])
}

mean(time_gower_mds)
mean(time_gower_candidate_1)
