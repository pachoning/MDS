source("tools/load_libraries.R")
source("tools/compute_accuracy.R")
source("tools/fast_MDS.R")


set.seed(12345)
n_obs = 10^3
x = data.frame(
  x1 = rnorm(n_obs),
  x2 = rnorm(n_obs),
  x3 = rnorm(n_obs),
  x4 = rnorm(n_obs),
  x5 = rnorm(n_obs),
  x6 = rnorm(n_obs)
)

results_fast_mds = fast_mds(
  x = x,
  number_coordinates = 2,
  metric = "euclidean",
  timeout = 1
)

head(x)
head(results_fast_mds)

results_compare_methods = compare_methods(
  mds_new_approach = results_fast_mds,
  x = x,
  metric = "euclidean",
  number_coordinates = 2
)

head(results_fast_mds, 8)
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

