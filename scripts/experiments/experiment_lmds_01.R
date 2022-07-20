#setwd("~/MEGA_tesis_cristian/mds_for_big_data/techincal/MDS")
source("tools/simulator.R")
source("final_methods/gower_interpolation_mds.R")
source("final_methods/emmanuel_mds.R")
source("final_methods/my_lmds.R")


results_folder_name <- "experiment_emmanuel_15"
experiment_label = "lmds_emmanuel_15"

sample_size = c(2*10^3, 5*10^3,
                10^4, 2*10^4, 5*10^4,
                10^5, 2*10^5, 5*10^5,
                10^6)
n_cols = c(10)

distribution_parameters <- list(list(var = rep(15, 5)))

scenarios <- list(sample_size = sample_size, n_cols = n_cols, distribution_parameters = distribution_parameters)
algorithms <- list(gower = gower_interpolation_mds, emmanuel = emmanuel_mds, lmds = my_lmds)

path <- file.path(getwd(), "data", "experiments", results_folder_name)
n_simulations <- 10
overwrite_simulations <- TRUE
n_sampling_points <- NA
l_lmds <- 1000
l_emmanuel <- 1000
l_gower <- 1000
n_cores <- 5 #5
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
  l_lmds = l_lmds,
  n_cores = n_cores,
  num_mds_dimesions = num_mds_dimesions,
  verbose = verbose
)
