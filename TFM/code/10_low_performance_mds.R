source("tools/load_libraries.R")

seq_wanted_mds = seq(from = 100, to = 3000, by = 100)
seq_wanted_distance = seq(from = 1000, to = 10000, by = 500)
n_cols = 5

# MDS
df_results = data.frame(
  n = rep(NA, length(seq_wanted_mds)),
  size_distance_matrix = rep(NA, length(seq_wanted_mds)),
  elpased_time_distance_matrix = rep(NA, length(seq_wanted_mds)),
  elapsed_time_mds = rep(NA, length(seq_wanted_mds))
)

i = 1
for(n_size in seq_wanted_mds){
  x <- matrix(rnorm(n_cols*n_size), ncol = n_cols)
  message(paste0("Doing for size: ", nrow(x)))
  # Time to get the distance matrix
  init_time_distance_matrix = proc.time()
  dist = cluster::daisy(
    x = x,
    metric = "euclidean"
  )
  elapsed_time_distance_matrix = proc.time() - init_time_distance_matrix
  elapsed_time_distance_matrix = as.numeric(elapsed_time_distance_matrix[3])
  
  
  
  # Size of the matrix distance
  distance_matrix_size = object.size(dist)/1024^2	
  
  
  # Time to get the MDS
  init_time_mds = proc.time()
  mds <- stats::cmdscale(
    d = dist,
    k = n_cols,
    eig = TRUE 
  )
  
  elapsed_time_mds = proc.time() - init_time_mds
  elapsed_time_mds = as.numeric(elapsed_time_mds[3])
  
  # Store
  df_results$n[i] = n_size
  df_results$size_distance_matrix[i] = distance_matrix_size
  df_results$elpased_time_distance_matrix[i] = elapsed_time_distance_matrix
  df_results$elapsed_time_mds[i] = elapsed_time_mds
  
  i = i + 1
}

View(df_results)

df_results %>% 
  ggplot(
    aes(
      x = n,
      y = elapsed_time_mds
    )
  ) +
  geom_line() + 
  labs(
    x ="sample size", 
    y = "Time (sec.)"
  )

# Memory
df_results_memory = data.frame(
  n = rep(NA, length(seq_wanted_distance)),
  size_distance_matrix = rep(NA, length(seq_wanted_distance)),
  elpased_time_distance_matrix = rep(NA, length(seq_wanted_distance))
)

i = 1
for(n_size in seq_wanted_distance){
  x <- matrix(rnorm(n_cols*n_size), ncol = n_cols)
  message(paste0("Doing for size: ", nrow(x)))
  # Time to get the distance matrix
  init_time_distance_matrix = proc.time()
  dist = cluster::daisy(
    x = x,
    metric = "euclidean"
  )
  elapsed_time_distance_matrix = proc.time() - init_time_distance_matrix
  elapsed_time_distance_matrix = as.numeric(elapsed_time_distance_matrix[3])
  
  
  
  # Size of the matrix distance
  distance_matrix_size = object.size(dist)/1024^2	

  
  # Store
  df_results_memory$n[i] = n_size
  df_results_memory$size_distance_matrix[i] = distance_matrix_size
  df_results_memory$elpased_time_distance_matrix[i] = elapsed_time_distance_matrix
  
  i = i + 1
}

df_results %>% 
  ggplot(
    aes(
      x = n,
      y = elapsed_time_mds
    )
  ) +
  geom_line() + 
  labs(
    x ="sample size", 
    y = "Time (sec.)"
  )

df_results_memory %>% 
  ggplot(
    aes(
      x = n,
      y = elpased_time_distance_matrix
    )
  ) +
  geom_line() + 
  labs(
    x ="sample size", 
    y = "Time (sec.)"
  )


df_results_memory %>% 
  ggplot(
    aes(
      x = n,
      y = size_distance_matrix
    )
  ) +
  geom_line() + 
  labs(
    x ="sample size", 
    y = "Memory consumed (MB)"
  )
