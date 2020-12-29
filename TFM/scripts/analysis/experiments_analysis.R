library(tidyverse)

data_path = file.path(getwd(), 'data')
load(file.path(data_path, "df_scenarios_full.RData"))
load(file.path(data_path, "df_time_full.RData"))
load(file.path(data_path, "df_correlation_full.RData"))


df_join_scenarios_time = df_time_full %>% 
  left_join(
    df_scenarios_full %>% select(-distribution_parameters, -mu, -sd, -processed_at, -computer_id),
    by = c("scenario_id" = "id")
  ) %>% 
  mutate(n_main_dimensions=as.factor(n_main_dimensions), 
         sample_size=as.factor(sample_size),
         n_cols = as.factor(n_cols))

df_join_scenarios_time %>% 
  group_by(method_name, n_main_dimensions, sample_size) %>% 
  summarise(mean_time = median(elapsed_time))%>%
  ggplot(aes(x=n_main_dimensions, y=mean_time, color=method_name, group=method_name)) +
  geom_line() +
  facet_wrap(. ~ sample_size, scales = "free")
 

df_correlation_full %>% View
