source("tools/simulator.R")
library(tidyverse)

sample_size = c(1000)
n_cols = c(10)

distribution_parameters = list(list(sd=c(NA)), list(sd=c(5)), list(sd=c(5,10)), list(sd=c(5,10, 15)), 
                               list(sd=c(5,10, 15, 20)), list(sd=c(5,10, 15, 25, 30)))
scenarios = list(sample_size=sample_size, n_cols=n_cols, distribution_parameters=distribution_parameters)


get_simulations(
  scenarios=scenarios,
  path=file.path(getwd(), 'data'),
  mds_methods = c(divide_conquer_mds, fast_mds, gower_interpolation_mds),
  n_simulations = 100,
  overwrite_simulations = TRUE,
  n_sampling_points = NA,
  largest_matrix_efficient_mds = 100,
  num_mds_dimesions = NA,
  verbose = TRUE
)


df_join = df_time %>% 
  left_join(
    df_scenarios %>% select(id, n_main_dimensions),
    by = c("scenario_id" = "id")
  )

df_join %>% 
  mutate(n_main_dimensions=as.factor(n_main_dimensions)) %>% 
  group_by(n_main_dimensions, method_name) %>% 
  summarise(mean_time = mean(elapsed_time))%>% 
  ggplot(aes(x=n_main_dimensions, y=mean_time, color=method_name, group=method_name)) +
  geom_line()
