library(tidyverse)
library(ggplot2)

# Load data ----
data_path = file.path(getwd(), 'data')
load(file.path(data_path, "df_scenarios_full.RData"))
load(file.path(data_path, "df_correlation_full.RData"))

# Manipulate data ----
scenario_identifier = c("sample_size", "n_cols", "n_main_dimensions", "sd_main")

# Avoid using scenarions which sample size is 10^6
df_scenarios_full_filtered = df_scenarios_full %>% 
  filter(sample_size != 10^6)

# Join scenarios and correlation
df_join_scenarios_correlation = df_scenarios_full_filtered %>% 
  select_at(c("id", scenario_identifier))   %>% 
  left_join(
    df_correlation_full %>% select(-n_main_dimensions),
    by = c("id" = "scenario_id")
  ) %>% 
  mutate(corr_1 = map2_dbl(.x = correlation_vector, .y = n_main_dimensions, .f = ~ifelse(.y >=0, .x[1], NA_real_)),
         corr_2 = map2_dbl(.x = correlation_vector, .y = n_main_dimensions, .f = ~ifelse(.y >=2, .x[2], NA_real_)),
         corr_3 = map2_dbl(.x = correlation_vector, .y = n_main_dimensions, .f = ~ifelse(.y >=3, .x[3], NA_real_)),
         corr_4 = map2_dbl(.x = correlation_vector, .y = n_main_dimensions, .f = ~ifelse(.y >=4, .x[4], NA_real_)))

# Analyse data  ----
# Plots
scenarios_with_main_dimesions = df_scenarios_full_filtered %>% 
  filter(n_main_dimensions > 0)

total_scenarios = nrow(scenarios_with_main_dimesions)
i = 1

while(i <= total_scenarios){
  current_scenario = scenarios_with_main_dimesions[i, ]
  scenario_id = current_scenario$id
  n_main_dimensions = current_scenario$n_main_dimensions
  
  correlation_data = df_join_scenarios_correlation[df_join_scenarios_correlation$id == scenario_id, ]
  cols_to_summarise = paste0("corr_", 1:n_main_dimensions)
  summarised_data = correlation_data %>% 
    group_by_at(c(scenario_identifier, "method_name")) %>% 
    summarise_at(.vars = cols_to_summarise, .funs = mean) %>% 
    pivot_longer(cols = all_of(cols_to_summarise))
  
  name_title = paste0("sample_size: ", current_scenario$sample_size, "; n_cols: ", current_scenario$n_cols, 
                      "; n_main_dim: ", current_scenario$n_main_dimensions, "; sd_main: ", current_scenario$sd_main)
  
  p = summarised_data %>% 
    ggplot(aes(x = name, y = value, color = method_name)) +
    geom_point() +
    ggtitle(name_title) +
    theme(plot.title = element_text(hjust = 0.5))
  
  png(file.path("plots", "correlation", paste0(name_title, ".png")))
  print(p)
  dev.off()
  
  i = i + 1
}

# Anova analysis
long_cols = paste0("corr_", 1:4)
data_wider = df_join_scenarios_correlation %>% 
  pivot_longer(all_of(long_cols), names_to = "dimension", values_to = "correlation") %>% 
  filter(!is.na(correlation)) %>% 
  mutate(
    n_cols = as.factor(n_cols),
    sample_size = as.factor(sample_size),
    n_main_dimensions = as.factor(n_main_dimensions),
    sd_main = as.factor(sd_main),
    method_name = as.factor(method_name),
    dimension = as.factor(dimension),
    log_correlation = log(correlation)
  )


# Measuring just the effect of the method
linear_model_method = lm(correlation ~ method_name, data = data_wider)
anova(linear_model_method)
summary(linear_model_method)

# Full Anova
linear_model = lm(correlation ~ n_cols + sample_size + n_main_dimensions + sd_main + method_name + dimension, data = data_wider)
anova(linear_model)
summary(linear_model)

# Selection of just one main dimension
data_wider_first_dimension = data_wider %>% filter(n_main_dimensions == 1)
linear_model_no_sd = lm(correlation ~ n_cols + sample_size + method_name, data = data_wider)
anova(linear_model_no_sd)
summary(linear_model_no_sd)
