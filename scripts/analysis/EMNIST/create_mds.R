source("final_methods/classical_mds.R")
source("final_methods/divide_conquer_mds.R")
source("final_methods/fast_mds.R")
source("final_methods/gower_interpolation_mds.R")


# Load data ----
load("data/EMNIST/all_data.RData")

# Manipulate data ----
all_data_pixels <- all_data_pixels/255

# Sample ----
sample_size <- NULL
labels_to_filter <- NULL
idx_target <- 1:length(target)

if (is.null(labels_to_filter)) {
  idx_target_specific <- idx_target
} else {
  idx_target_specific <- idx_target[target %in% labels_to_filter]
}

if (is.null(sample_size)) {
  sample_size <- length(idx_target_specific)
}

idx_sample <- sample(x = idx_target_specific, size = pmin(sample_size, length(idx_target_specific)))
idx_sample <- sort(idx_sample)
target <- target[idx_sample]
all_data_pixels <- all_data_pixels[idx_sample, ]

# Fast ----
init_proc <- proc.time()
fast_results <- fast_mds(x = all_data_pixels, l = 200, s = 20, k = 2)
elapsed_time_fast <- proc.time() - init_proc
elapsed_time_fast <- elapsed_time_fast[3]
mds_fast <- fast_results$points
df_mds_fast <- data.frame(x =  mds_fast[, 1], y = mds_fast[, 2], target = target)
eigen_fast <- fast_results$eigen
GOF_fast <- fast_results$GOF
save(df_mds_fast, eigen_fast, GOF_fast, elapsed_time_fast, file = "data/EMNIST/fast.RData")

# Gower ----
init_proc <- proc.time()
gower_results <- gower_interpolation_mds(x = all_data_pixels, l = 200, k = 2)
elapsed_time_gower <- proc.time() - init_proc
elapsed_time_gower <- elapsed_time_gower[3]
mds_gower <- gower_results$points
df_mds_gower <- data.frame(x =  mds_gower[, 1], y = mds_gower[, 2], target = target)
eigen_gower <- gower_results$eigen
GOF_gower <- gower_results$GOF
save(df_mds_gower, eigen_gower, GOF_gower, elapsed_time_gower, file = "data/EMNIST/gower.RData")

# Divide and conquer ----
init_proc <- proc.time()
divide_results <- divide_conquer_mds(x = all_data_pixels, l = 200, tie = 20, k = 2)
elapsed_time_divide <- proc.time() - init_proc
elapsed_time_divide <- elapsed_time_divide[3]
mds_divide <- divide_results$points
df_mds_divide <- data.frame(x =  mds_divide[, 1], y = mds_divide[, 2], target = target)
eigen_divide <- divide_results$eigen
GOF_divide <- divide_results$GOF
save(df_mds_divide, eigen_divide, GOF_divide, elapsed_time_divide, file = "data/EMNIST/divide.RData")
