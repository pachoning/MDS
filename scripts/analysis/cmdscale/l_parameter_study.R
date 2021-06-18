dominant_dimension <- rep(15, 5)
n_dominant_dimension <- length(dominant_dimension)
rest_dimensions <- rep(1, 5)
var <- c(dominant_dimension, rest_dimensions)
l_values <- seq(from = 100, to = 5000, by = 100)

df_cost_time <- data.frame(
  l = numeric(0),
  event = character(0),
  time = numeric(0))

i_iter <- 1
for (l in l_values) {

  msg <- paste0("Working on l: ", l)
  message(msg)
  
  sample_size <- 10^5
  x <- mapply(rnorm, sd = sqrt(var), MoreArgs = list(n = l))
  init_time <- proc.time()
  dist_matrix <- dist(x)
  end_time <- proc.time() - init_time 
  elapsed_time_dist <- proc.time() - init_time 
  elapsed_time_dist <- as.numeric(elapsed_time_dist[3])
  
  init_time <- proc.time()
  results <- cmdscale(d = dist_matrix, k = n_dominant_dimension, eig = TRUE)
  end_time <- proc.time() - init_time
  elapsed_time_cmdscale <- proc.time() - init_time
  elapsed_time_cmdscale <- as.numeric(elapsed_time_cmdscale[3])
  
  current_iter_results <- data.frame(l = c(l, l),
                                     event = c("dist_matrix", "cmdscale"),
                                     time = c(elapsed_time_dist, elapsed_time_cmdscale))
  df_cost_time <- rbind(df_cost_time, current_iter_results)
  
  save(df_cost_time, file = file.path("data", "cmdscale", "df_cost_time.RData"))
  i_iter <- i_iter + 1
  
}
