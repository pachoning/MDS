source("final_methods/divide_conquer_mds.R")
source("final_methods/fast_mds.R")
source("final_methods/gower_interpolation_mds.R")

# Load data ----
load("data/EMNIST/all_data.RData")

# Manipulate data ----
all_data_pixels <- all_data_pixels/255
n_cores <- 5
n_main_dim <- 2

is_80_reached <- FALSE
i_iter <- 1

df_emnist_gof <- data.frame(
  method_name = character(0),
  k = numeric(0),
  GOF = numeric(0)
)

while (!is_80_reached) {
  
  message(paste0("Starting iteration: ", i_iter, ". Num. main dimensions: ", n_main_dim))
  
  # Divide
  message(paste0("\tStarting Divide"))
  divide_results <- divide_conquer_mds(x = all_data_pixels, l = 400, tie = 100, k = n_main_dim, n_cores = n_cores)
  divide_GOF <- divide_results$GOF[1]
  
  # Gower
  message(paste0("\tStarting Gower"))
  gower_results <- gower_interpolation_mds(x = all_data_pixels, l = 200, k = n_main_dim, n_cores = n_cores)
  gower_GOF <- gower_results$GOF[1]
  
  # Fast
  message(paste0("\tStarting Fast"))
  fast_results <- fast_mds(x = all_data_pixels, l = 1000, s = 100, k = n_main_dim, n_cores = n_cores)
  fast_GOF <- fast_results$GOF[1]
  
  df_temp_gof <- data.frame(
    method_name = c("divide", "gower", "fast"),
    k = c(n_main_dim, n_main_dim, n_main_dim),
    GOF = c(divide_GOF, gower_GOF, fast_GOF)
  )
  df_emnist_gof <- rbind(df_emnist_gof, df_temp_gof)
  
  i_iter <- i_iter + 1
  n_main_dim <- n_main_dim + 1
  save(df_emnist_gof, file = "data/EMNIST/df_emnist_gof.RData")
  
  if (min(c(divide_GOF, gower_GOF, fast_GOF)) > 0.8) {
    is_80_reached <- TRUE
  }
}
