library(tidyverse)
library(ggplot2)

# Load data ----
data_path = file.path(getwd(), 'data')
load(file.path(data_path, "df_scenarios_full.RData"))
load(file.path(data_path, "df_correlation_full.RData"))

# Manipulate data ----
scenario_identifier = c("sample_size", "n_cols", "n_main_dimensions", "sd_main")

# Avoid using scenarions which sample size is 10^6
df_scenarios_full_filtered = df_scenarios_full

# Join scenarios and correlation
df_join_scenarios_correlation = df_scenarios_full_filtered %>% 
  select_at(c("id", scenario_identifier))   %>% 
  left_join(
    df_correlation_full %>% select(-n_main_dimensions),
    by = c("id" = "scenario_id")
  ) %>% 
  mutate(dim_1 = map2_dbl(.x = correlation_vector, .y = n_main_dimensions, .f = ~ifelse(.y >=0, .x[1], NA_real_)),
         dim_2 = map2_dbl(.x = correlation_vector, .y = n_main_dimensions, .f = ~ifelse(.y >=2, .x[2], NA_real_)),
         dim_3 = map2_dbl(.x = correlation_vector, .y = n_main_dimensions, .f = ~ifelse(.y >=3, .x[3], NA_real_)),
         dim_4 = map2_dbl(.x = correlation_vector, .y = n_main_dimensions, .f = ~ifelse(.y >=4, .x[4], NA_real_))) %>% 
  select(-correlation_vector) %>% 
  pivot_longer(cols = starts_with("dim"), names_to = "dim", values_to = "correlation") %>% 
  filter(!is.na(correlation)) %>% 
  mutate(log_1_correlation = log(1-correlation))

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
  
  name_title = paste0("sample_size: ", current_scenario$sample_size, "; n_cols: ", current_scenario$n_cols, 
                      "; n_main_dim: ", current_scenario$n_main_dimensions, "; sd_main: ", current_scenario$sd_main)
  
  p = correlation_data %>% 
    ggplot(aes(x = method_name, y = correlation)) + 
    geom_boxplot() + 
    facet_wrap(. ~ dim, ncol = n_main_dimensions) +
    ggtitle(name_title) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  png(file.path("plots", "correlation", "plot_correlation", paste0(name_title, ".png")))
  print(p)
  dev.off()
  
  p = correlation_data %>% 
    ggplot(aes(x = method_name, y = log_1_correlation, fill = method_name)) + 
    geom_boxplot() + 
    facet_wrap(. ~ dim, ncol = n_main_dimensions) +
    ggtitle(name_title) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none")
  
  png(file.path("plots", "correlation", "plot_log_correlation", paste0(name_title, ".png")))
  print(p)
  dev.off()
  
  i = i + 1
}

# Anova analysis
# Measuring just the effect of the method
linear_model_method = lm(correlation ~ method_name, data = df_join_scenarios_correlation)
anova(linear_model_method)
summary(linear_model_method)

# Add sample size and interaction between method_name and sample size
linear_model_method_sample_size = df_join_scenarios_correlation %>% mutate(sample_size = as.factor(sample_size))%>% 
  lm(correlation ~ method_name + sample_size + method_name:sample_size, data = .)
anova(linear_model_method_sample_size)
summary(linear_model_method_sample_size)
