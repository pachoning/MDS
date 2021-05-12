n_col <- 5
v_sample_size <- seq(from = 3*n_col, to = 10000, by = 100)
mean_val <- rep(0, n_col)
sd_val <- rep(1, n_col)

df_info_process <- data.frame(
  sample_size = v_sample_size, 
  time_distance = NA, 
  time_mds = NA,
  memory_distance = NA,
  memory_mds = NA)

i = 1
for(sample_size in v_sample_size) {
  message(paste0("Processing ", i, " out of ", length(v_sample_size)))
  X <- mapply(n = sample_size, rnorm, mean = mean_val, sd = sd_val)
  
  init_time <- proc.time()
  dist_X <- dist(X)
  elapsed_time <- proc.time() - init_time
  size_dist_X <- object.size(dist_X)
  df_info_process[i, 'time_distance'] <- elapsed_time[3]
  df_info_process[i, 'memory_distance'] <- size_dist_X
  
  init_time <- proc.time()
  mds_result <- cmdscale(d = dist_X, k = n_col)
  elapsed_time <- proc.time() - init_time
  size_dist_X <- object.size(mds_result)
  df_info_process[i, 'time_mds'] <- elapsed_time[3]
  df_info_process[i, 'memory_mds'] <- size_dist_X
  
  i = i + 1
}

df_info_process_filter <- df_info_process[!is.na(df_info_process$memory_distance), ]
View(df_info_process_filter)
save(df_info_process_filter, file="data/other/df_info_process_filter.RData")
