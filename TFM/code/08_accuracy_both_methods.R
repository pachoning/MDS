# ftp://ftp.auckland.ac.nz/pub/software/CRAN/doc/packages/micEcdat.pdf
library(Ecdat)
source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")
source("tools/fast_mds.R")
source("tools/compute_accuracy.R")


load("data/Bike-Sharing-Dataset/df_split.RData")



set.seed(12345)

# 01 Random example ----
# Creation od data set
n_obs = 3*10^3
x = data.frame(
  x1 = rnorm(n_obs),
  x2 = rnorm(n_obs),
  x3 = rnorm(n_obs),
  x4 = rnorm(n_obs),
  x5 = rnorm(n_obs),
  x6 = rnorm(n_obs)
)


# Divide and conquer MDS
metric = "euclidean"
mds_divide_conquer = divide_conquer_mds(
  x = x,
  groups =  sample(x = 3, size = nrow(x), replace = TRUE),
  # groups = c(rep(1,25), 26),
  number_coordinates = 2,
  metric = metric
)


results_compare_divide_conquer = compare_methods(
  mds_new_approach = mds_divide_conquer$mds,
  x = x,
  metric = metric,
  number_coordinates = 2
)
head(mds_divide_conquer$mds, 8)
head(results_compare_divide_conquer$mds_classical_transformed, 8)

# Plot coordinates
results_compare_divide_conquer$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_divide_conquer$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_divide_conquer$distance_between_coordinates)

################################################################################

# Fast MDS
mds_fast = fast_mds(
  x = x,
  number_coordinates = 2,
  metric = metric,
  timeout = 1
)

results_compare_fast = compare_methods(
  mds_new_approach = mds_fast,
  x = x,
  metric = metric,
  number_coordinates = 2
)


head(mds_fast, 8)
head(results_compare_fast$mds_classical_transformed, 8)


# Plot coordinates
results_compare_fast$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_fast$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_fast$distance_between_coordinates)



# 02 iris ----
metric = "euclidean"
x = iris[, -5]

mds_divide_conquer = divide_conquer_mds(
  x = x,
  groups =  sample(x = 2, size = nrow(x), replace = TRUE),
  # groups = c(rep(1,25), 26),
  number_coordinates = 2,
  metric = metric
)



results_compare_divide_conquer = compare_methods(
  mds_new_approach = mds_divide_conquer$mds,
  x = x,
  metric = "euclidean",
  number_coordinates = 2
)

head(mds_divide_conquer$mds, 8)
head(results_compare_divide_conquer$mds_classical_transformed, 8)


# Plot coordinates
results_compare_divide_conquer$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_divide_conquer$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_divide_conquer$distance_between_coordinates)




# 03 Bike sharing data set ----
metric = "gower"
x = df_split %>%
  slice(1:4000) %>%
  select(
    -instant,
    -dteday,
    -group_member
  )



# Divide and conquer MDS
mds_divide_conquer = divide_conquer_mds(
  x = x,
  groups =  sample(x = 5, size = nrow(x), replace = TRUE),
  number_coordinates = 2,
  metric = metric
)


results_compare_divide_conquer = compare_methods(
  mds_new_approach = mds_divide_conquer$mds,
  x = x,
  metric = metric,
  number_coordinates = 2
)


head(mds_divide_conquer$mds, 8)
head(results_compare_divide_conquer$mds_classical_transformed, 8)

# Plot coordinates
results_compare_divide_conquer$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_divide_conquer$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_divide_conquer$distance_between_coordinates)


################################################################################

# Fast MDS
mds_fast = fast_mds(
  x = x,
  number_coordinates = 2,
  metric = metric,
  timeout = 1
)

results_compare_fast = compare_methods(
  mds_new_approach = mds_fast,
  x = x,
  metric = "gower",
  number_coordinates = 2
)


head(mds_fast, 8)
head(results_compare_fast$mds_classical_transformed, 8)


# Plot coordinates
results_compare_fast$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_fast$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_fast$distance_between_coordinates)



# 04 Budget Food ----

set.seed(12345)
data(BudgetFood)
nrow(BudgetFood)
ind = sample(x = nrow(BudgetFood), size = 3000, replace = FALSE)
x = BudgetFood %>% slice(ind) %>% select(-sex, -town) 
metric = "euclidean"


# Divide and conquer MDS
mds_divide_conquer = divide_conquer_mds(
  x = x,
  groups =  sample(x = 5, size = nrow(x), replace = TRUE),
  number_coordinates = 2,
  metric = metric
)


results_compare_divide_conquer = compare_methods(
  mds_new_approach = mds_divide_conquer$mds,
  x = x,
  metric = metric,
  number_coordinates = 2
)


head(mds_divide_conquer$mds, 8)
head(results_compare_divide_conquer$mds_classical_transformed, 8)

# Plot coordinates
results_compare_divide_conquer$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_divide_conquer$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_divide_conquer$distance_between_coordinates)



rdist(
  results_compare_divide_conquer$mds_classical_transformed,
  procrustes_result$X.new
)
