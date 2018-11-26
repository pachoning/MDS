source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")
source("tools/compute_accuracy.R")

set.seed(12345)

# 01 Very random example ----
n_obs = 50
x = data.frame(
  x1 = rnorm(n_obs),
  x2 = rnorm(n_obs),
  x3 = rnorm(n_obs),
  x4 = rnorm(n_obs),
  x5 = rnorm(n_obs),
  x6 = rnorm(n_obs)
)


# Divide and conquer MDS
mds_divide_conquer = divide_conquer_mds(
  x = x,
  groups =  sample(x = 2, size = nrow(x), replace = TRUE),
  # groups = c(rep(1,25), 26),
  number_coordinates = 2,
  metric = "euclidean"
)


results_compare_methods = compare_methods(
  mds_new_approach = mds_divide_conquer$mds,
  x = x,
  metric = "euclidean",
  number_coordinates = 2
)
head(mds_divide_conquer$mds, 8)
head(results_compare_methods$mds_classical_transformed, 8)
head(results_compare_methods$mds_classical, 8)

# Plot coordinates
results_compare_methods$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_methods$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_methods$distance_between_coordinates)



# 02 Iris data set ----
x = iris[, -5]


mds_divide_conquer = divide_conquer_mds(
  x = x,
  groups =  sample(x = 2, size = nrow(x), replace = TRUE),
  # groups = c(rep(1,25), 26),
  number_coordinates = 2,
  metric = "euclidean"
)



results_compare_methods = compare_methods(
  mds_new_approach = mds_divide_conquer$mds,
  x = x,
  metric = "euclidean",
  number_coordinates = 2
)

head(mds_divide_conquer$mds, 8)
head(results_compare_methods$mds_classical_transformed, 8)
head(results_compare_methods$mds_classical, 8)

# Plot coordinates
results_compare_methods$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_methods$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_methods$distance_between_coordinates)


# 03 MNIST data set ----



