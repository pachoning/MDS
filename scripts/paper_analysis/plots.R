library(tidyverse)
library(ggplot2)

# Load data ----
data_path <- file.path(getwd(), 'data')
load(file.path(data_path, "df_scenarios_full.RData"))
load(file.path(data_path, "df_correlation_full.RData"))


# Manipulate data ----
scenario_identifier <- c("sample_size", "n_cols", "n_main_dimensions", "sd_main")

# Join scenarios and correlation
df_join_scenarios_correlation <- scenarios_with_main_dimesions %>% 
  select_at(c("id", scenario_identifier))   %>% 
  left_join(
    df_correlation_full %>% select(-n_main_dimensions),
    by = c("id" = "scenario_id")
  ) %>% unnest(correlation_vector) %>% 
  rename(correlation = correlation_vector) %>% 
  group_by(id, method_name, num_sim) %>% 
  mutate(
    dim = paste0("dim_", 1:n()),
  )


# Results/plots ----
filter <- df_join_scenarios_correlation$sample_size == 1000 & 
  df_join_scenarios_correlation$n_main_dimensions == 3 &
  df_join_scenarios_correlation$method_name == "fast"

scenario_selection <- df_join_scenarios_correlation[filter, ]
