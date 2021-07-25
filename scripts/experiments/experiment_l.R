source("tools/load_libraries.R")
source("tools/simulator.R")
source("final_methods/divide_conquer_mds.R")
source("final_methods/gower_interpolation_mds.R")
source("final_methods/fast_mds.R")


results_folder_name <- "experiment_110"
experiment_label <- "l_experiment"

sample_size <- c(10^6)
n_cols <- c(100)

distribution_parameters <- list(list(var = rep(15, 10)))

scenarios <- list(sample_size = sample_size, 
                 n_cols = n_cols, 
                 distribution_parameters = distribution_parameters)

algorithms <- list(divide = divide_conquer_mds,
                  gower = gower_interpolation_mds,
                  fast = fast_mds)

path <- file.path(getwd(), "data", "experiments", results_folder_name)
n_simulations <- 1
overwrite_simulations <- FALSE
n_sampling_points <- NA
l_default <- 1500
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
  l_divide = l_default,
  l_gower = l_default,
  l_fast = l_default,
  n_cores = n_cores,
  num_mds_dimesions = num_mds_dimesions,
  verbose = verbose
)
