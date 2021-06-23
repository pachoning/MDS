source("tools/simulator.R")

results_folder_name = "experiment_02"
experiment_label = "10^4"

sample_size = c(10000)
n_cols = c(10, 100)

distribution_parameters = list(list(var=15), list(var=rep(15, 2)), list(var=rep(15, 3)), list(var=rep(15, 4)),
                               list(var=rep(15, 5)), list(var=rep(15, 6)), list(var=rep(15, 7)), 
                               list(var=rep(15, 8)), list(var=rep(15, 9)), list(var=rep(15, 10)))

scenarios = list(sample_size=sample_size, n_cols=n_cols, distribution_parameters=distribution_parameters)

get_simulations(
  experiment_label=experiment_label, 
  scenarios=scenarios,
  path=file.path(getwd(), "data", "experiments", results_folder_name),
  mds_methods_names=c("divide_conquer", "fast", "gower"),
  n_simulations=100,
  overwrite_simulations=FALSE,
  n_sampling_points=NA,
  largest_matrix_efficient_mds=200,
  num_mds_dimesions=NA,
  verbose=TRUE
)
