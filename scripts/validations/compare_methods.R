library(tidyverse)
source("tools/load_libraries.R")
source("tools/procrustes.R")
source("final_methods/divide_conquer_stepwise_mds.R")
source("tools/gower_by_row.R")
source("final_methods/gower_interpolation_mds.R")


sample_size <- c(10^3, 10^4, 10^5, 10^6)
n_sim <- 20
methods <- c("gower_by_row", "gower_interpolation_mds")
total_length <- length(sample_size) * n_sim * length(methods)
l <- 200

info <- data.frame(method = rep(NA, total_length),
                 sample_size = rep(NA, total_length),
                 elapsed_time = rep(NA, total_length))
iter <- 1
for (ss in sample_size) {
  for (i in 1:n_sim) {
    msg <- paste0("Iteration: ", iter," out of: ", total_length)
    message(msg)
    sd_vector <- rep(10, 10)
    mean_vector <- rep(0, 10)
    n_cols <- length(sd_vector)
    x <- mapply(rnorm, n = ss, mean_vector) %*% diag(sd_vector)
    
    for (name in methods) {
      if (name == "gower_by_row") {
        method = gower_by_row
      } else {
        method = gower_interpolation_mds
      }
      init_time <- proc.time()
      results <- method(x = x, l = l, k = n_cols)
      end_time <- proc.time()
      elapsed <- (end_time - init_time)[3]
      info$method[iter] <- name
      info$sample_size[iter] <- ss
      info$elapsed_time[iter] <- elapsed
      iter <- iter + 1
    }
  }
}


info %>% 
  filter(!is.na(sample_size)) %>% 
  group_by(sample_size, method) %>% 
  summarise(mean(elapsed_time)) %>% View
