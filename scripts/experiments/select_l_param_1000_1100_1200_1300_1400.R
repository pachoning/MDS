source("final_methods/fast_mds.R")
source("final_methods/interpolation_mds.R")
source("final_methods/divide_conquer_mds.R")
source("final_methods/landmark_mds.R")
source("final_methods/pivot_mds.R")
source("final_methods/reduced_mds.R")
source("tools/procrustes.R")
source("tools/new_simulator.R")


simulator(
  methods = list(
    #divide_conquer_mds = divide_conquer_mds,
    #interpolation_mds = interpolation_mds,
    #fast_mds = fast_mds,
    #landmark_mds = landmark_mds
    #pivot_mds = pivot_mds
    reduced_mds = reduced_mds
    
  ),
  n_samples = list(1000000),
  variances = list(
    c(rep.int(15, 10), rep.int(1, 90))
  ),
  l_global = as.list(c(1000, 1100, 1200, 1300, 1400)),
  n_sim = 20,
  path = file.path("data", "l_param_reduced_1000_1100_1200_1300_1400")
)
