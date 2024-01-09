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
    divide_conquer_mds = divide_conquer_mds,
    interpolation_mds = interpolation_mds,
    fast_mds = fast_mds,
    landmark_mds = landmark_mds,
    pivot_mds = pivot_mds,
    reduced_mds = reduced_mds
    
  ),
  n_samples = list(5000, 10000, 20000, 100000, 1000000),
  variances = list(
    c(rep.int(15, 2), rep.int(1, 8)),
    rep.int(15, 10),
    c(rep.int(15, 2), rep.int(1, 98)),
    c(rep.int(15, 10), rep.int(1, 90))
  ),
  l_methods =list(
    divide_conquer_mds = 150,
    interpolation_mds = 150,
    fast_mds = 400,
    landmark_mds = 150,
    pivot_mds = 150,
    reduced_mds = 50
  ),
  n_sim = 20,
  path = file.path("data", "full_experiments")
)
