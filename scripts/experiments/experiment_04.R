source("tools/simulator.R")
source("final_methods/divide_conquer_mds.R")
source("final_methods/gower_interpolation_mds.R")
source("final_methods/fast_mds.R")

results_folder_name <- "experiment_04"
experiment_label <- "experiment_10^5"

sample_size <- c(10^5)
n_cols <- c(10, 100)

distribution_parameters <- list(list(var = 15), list(var = rep(15, 2)), list(var = rep(15, 3)), 
                                list(var = rep(15, 4)), list(var = rep(15, 5)), list(var = rep(15, 6)), 
                                list(var = rep(15, 7)), list(var = rep(15, 8)), list(var = rep(15, 9)), 
                                list(var = rep(15, 10)))

scenarios <- list(sample_size = sample_size, n_cols = n_cols, distribution_parameters = distribution_parameters)
algorithms <- list(divide = divide_conquer_mds, gower = gower_interpolation_mds, fast = fast_mds)

path <- file.path(getwd(), "data", "experiments", results_folder_name)
n_simulations <- 100
overwrite_simulations <- FALSE
n_sampling_points <- NA
l_divide <- 500
l_gower <- 1500
l_fast <- 1500
n_cores <- 5
num_mds_dimesions <- NA
verbose <- TRUE

get_simulations(
  experiment_label = experiment_label, 
  scenarios = scenarios,
  path = path,
  algorithms = algorithms,
  n_simulations = n_simulations,
  overwrite_simulations = overwrite_simulations,
  n_sampling_points = n_sampling_points,
  l_divide = l_divide,
  l_gower = l_gower,
  l_fast = l_fast,
  n_cores = n_cores,
  num_mds_dimesions = num_mds_dimesions,
  verbose = verbose
)
