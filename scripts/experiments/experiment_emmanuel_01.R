source("tools/simulator.R")
source("final_methods/gower_interpolation_mds.R")
source("final_methods/emmanuel_mds.R")


results_folder_name <- "experiment_emmanuel_01"
experiment_label = "emmanuel_01"

sample_size = c(2*10^3, 5*10^3, 
                10^4, 2*10^4, 5*10^4, 
                10^5, 2*10^5, 5*10^4, 
                10^6)
n_cols = c(10)

distribution_parameters <- list(list(var = rep(15, 5)))

scenarios <- list(sample_size = sample_size, n_cols = n_cols, distribution_parameters = distribution_parameters)
algorithms <- list(emmanuel = emmanuel_mds, gower = gower_interpolation_mds)

path <- file.path(getwd(), "data", "experiments", results_folder_name)
n_simulations <- 1
overwrite_simulations <- FALSE
n_sampling_points <- NA
l_fast <- 1000
l_emmanuel <- 1000
l_gower <- 1000
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
  l_gower = l_gower,
  l_emmanuel = l_emmanuel,
  n_cores = n_cores,
  num_mds_dimesions = num_mds_dimesions,
  verbose = verbose
)
