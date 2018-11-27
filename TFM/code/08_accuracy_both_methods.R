# ftp://ftp.auckland.ac.nz/pub/software/CRAN/doc/packages/micEcdat.pdf
library(Ecdat)
source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")
source("tools/fast_mds.R")
source("tools/compute_accuracy.R")


load("data/Bike-Sharing-Dataset/df_split.RData")


load("/Users/Cristian/Documents/MESIO/TFM/check_point_2018_11_26_19_09_00.Rproj")
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

x <- as.data.frame(matrix(rnorm(5*n_obs),nrow=n_obs))


# Divide and conquer MDS
metric = "euclidean"
mds_divide_conquer_random = divide_conquer_mds(
  x = x,
  groups =  sample(x = 3, size = nrow(x), replace = TRUE),
  # groups = c(rep(1,25), 26),
  number_coordinates = 2,
  metric = metric
)

results_classical_mds_random = classical_mds(
  x = x,
  number_coordinates = 2,
  metric = metric
)

results_compare_divide_conquer_random = compare_methods(
  mds_new_approach = mds_divide_conquer_random$mds,
  mds_classical = results_classical_mds_random
)

head(mds_divide_conquer_random$mds, 8)
head(results_compare_divide_conquer_random$mds_classical_transformed, 8)

head(results_compare_divide_conquer_random$df_both_mds_labels)

# Plot coordinates
results_compare_divide_conquer_random$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_divide_conquer_random$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_divide_conquer_random$distance_between_coordinates)


df1 = results_compare_divide_conquer_random$df_both_mds_labels %>% 
  filter(type == "classical") %>% 
  select(V1) %>% 
  pull()

df1

df2 = results_compare_divide_conquer_random$df_both_mds_labels %>% 
  filter(type == "classical") %>% 
  select(V2)

df2


df1n = results_compare_divide_conquer_random$df_both_mds_labels %>% 
  filter(type == "new_approach") %>% 
  select(V1) %>% 

df1n 

df2n = results_compare_divide_conquer_random$df_both_mds_labels %>% 
  filter(type == "new_approach") %>% 
  select(V2)

df2n

plot(df1[[1]], df1n[[1]])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(df2[[1]], df2n[[1]])
abline(a = 0, b = 1, col = 2, lwd = 2)


################################################################################

# Fast MDS
mds_fast_random = fast_mds(
  mds_classical = x,
  number_coordinates = 2,
  metric = metric,
  timeout = 1
)

results_compare_fast_random = compare_methods(
  mds_new_approach = mds_fast_random,
  mds_classical = results_classical_mds_random
)


head(mds_fast_random, 8)
head(results_compare_fast_random$mds_classical_transformed, 8)


# Plot coordinates
results_compare_fast_random$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_fast_random$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_fast_random$distance_between_coordinates)



# 02 iris ----
metric = "euclidean"
x = iris[, -5]

mds_divide_conquer_iris = divide_conquer_mds(
  x = x,
  groups =  sample(x = 2, size = nrow(x), replace = TRUE),
  # groups = c(rep(1,25), 26),
  number_coordinates = 2,
  metric = metric
)


results_classical_mds_iris = classical_mds(
  x = x,
  number_coordinates = 2,
  metric = metric
)


results_compare_divide_conquer_iris = compare_methods(
  mds_new_approach = mds_divide_conquer_iris$mds,
  mds_classical = results_classical_mds_iris
)

head(mds_divide_conquer_iris$mds, 8)
head(results_compare_divide_conquer_iris$mds_classical_transformed, 8)


# Plot coordinates
results_compare_divide_conquer_iris$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_divide_conquer_iris$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_divide_conquer_iris$distance_between_coordinates)




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
mds_divide_conquer_bike = divide_conquer_mds(
  x = x,
  groups =  sample(x = 5, size = nrow(x), replace = TRUE),
  number_coordinates = 2,
  metric = metric
)


results_classical_mds_bike = classical_mds(
  x = x,
  number_coordinates = 2,
  metric = metric
)

results_compare_divide_conquer_bike = compare_methods(
  mds_new_approach = mds_divide_conquer_bike$mds,
  mds_classical = results_classical_mds_bike
)


head(mds_divide_conquer_bike$mds, 8)
head(results_compare_divide_conquer_bike$mds_classical_transformed, 8)

# Plot coordinates
results_compare_divide_conquer_bike$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_divide_conquer_bike$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_divide_conquer_bike$distance_between_coordinates)


################################################################################

# Fast MDS
mds_fast_bike = fast_mds(
  x = x,
  number_coordinates = 2,
  metric = metric,
  timeout = 1
)

results_compare_fast_bike = compare_methods(
  mds_new_approach = mds_fast_bike,
  mds_classical = mds_fast_bike
)


head(mds_fast_bike, 8)
head(results_compare_fast_bike$mds_classical_transformed, 8)


# Plot coordinates
results_compare_fast_bike$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_fast_bike$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_fast_bike$distance_between_coordinates)



# 04 Budget Food ----

set.seed(12345)
data(BudgetFood)
nrow(BudgetFood)
ind = sample(x = nrow(BudgetFood), size = 3000, replace = FALSE)
x = BudgetFood %>% slice(ind) %>% select(-sex, -town) 
metric = "euclidean"


# Divide and conquer MDS
mds_divide_conquer_budget = divide_conquer_mds(
  x = x,
  groups =  sample(x = 5, size = nrow(x), replace = TRUE),
  number_coordinates = 2,
  metric = metric
)

results_classical_mds_budget = classical_mds(
  x = x,
  number_coordinates = 2,
  metric = metric
)


results_compare_divide_conquer_budget = compare_methods(
  mds_new_approach = mds_divide_conquer_budget$mds,
  mds_classical = results_classical_mds_budget
)


head(mds_divide_conquer_budget$mds, 8)
head(results_compare_divide_conquer_budget$mds_classical_transformed, 8)

# Plot coordinates
results_compare_divide_conquer_budget$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_divide_conquer_budget$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_divide_conquer_budget$distance_between_coordinates)


################################################################################

# Fast MDS
mds_fast_budget = fast_mds(
  x = x,
  number_coordinates = 2,
  metric = metric,
  timeout = 1
)

results_compare_fast_budget = compare_methods(
  mds_new_approach = mds_fast_budget,
  mds_classical = mds_fast_budget
)


head(mds_fast_budget, 8)
head(results_compare_fast_budget$mds_classical_transformed, 8)


# Plot coordinates
results_compare_fast_budget$df_both_mds_labels[1:20, ] %>% 
  ggplot(aes(x = V1, y = V2, color = type, group = type)) + 
  geom_text(aes(label=label),hjust=0, vjust=0)


# Plot the error
ggplot(
  data.frame(
    error = results_compare_fast_budget$distance_between_coordinates
  ), 
  aes(error)
) +
  geom_density()

summary(results_compare_fast_bike$distance_between_coordinates)

save.image(file = "show_to_pedro.Rproj")

