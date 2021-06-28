source("tools/load_libraries.R")
source("tools/simulator.R")
source("tools/classical_mds.R")
source("final_methods/classical_mds.R")

source("final_methods/divide_conquer_mds.R")
source("final_methods/divide_conquer_small_mds.R")

source("final_methods/gower_interpolation_mds.R")
source("final_methods/gower_interpolation_rowwise_mds.R")
source("final_methods/gower_interpolation_kgropus_mds.R")
source("final_methods/gower_interpolation_dist_mds.R")
source("final_methods/gower_interpolation_big_mds.R")

source("final_methods/fast_mds.R")


results_folder_name = "experiment_99"
experiment_label = "using_do_call"

sample_size = c(10^4)
n_cols = c(100)

distribution_parameters = list(list(var=c(15, 15, 15, 15, 15)))

scenarios = list(sample_size = sample_size, 
                 n_cols = n_cols, 
                 distribution_parameters = distribution_parameters)

algorithms = list(#divide = divide_conquer_mds,
                  #divide_small = divide_conquer_small_mds,
                  
                  gower = gower_interpolation_mds,
                  gower_rowwise = gower_interpolation_rowwise_mds,
                  gower_kgroups = gower_interpolation_kgroups_mds,
                  gower_dist = gower_interpolation_dist_mds,
                  gower_big = gower_interpolation_big_mds,
                  
                  fast_mds = fast_mds)


path = file.path(getwd(), "data", "experiments", results_folder_name)
n_simulations = 3
overwrite_simulations = FALSE
n_sampling_points = NA
largest_matrix_efficient_mds = 1500
num_mds_dimesions = NA
verbose = TRUE

get_simulations(
  experiment_label = experiment_label, 
  scenarios = scenarios,
  path = path,
  algorithms = algorithms,
  n_simulations = n_simulations,
  overwrite_simulations = overwrite_simulations,
  n_sampling_points = n_sampling_points,
  largest_matrix_efficient_mds = largest_matrix_efficient_mds,
  num_mds_dimesions = num_mds_dimesions,
  verbose = verbose
)
