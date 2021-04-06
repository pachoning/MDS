library(tidyverse)
library(ggplot2)

# Load data ----
data_path = file.path(getwd(), 'data')
load(file.path(data_path, "df_scenarios_full.RData"))
load(file.path(data_path, "df_eigenvalues_full.RData"))

# Manipulate data ----
scenario_identifier = c("sample_size", "n_cols", "n_main_dimensions", "sd_main")

# Avoid using scenarions which sample size is 10^6
df_scenarios_full_filtered = df_scenarios_full

# Join scenarios and eigenvalues
scenarios_with_main_dimensions = df_scenarios_full_filtered %>% 
  filter(n_main_dimensions > 0)
total_scenarios = nrow(scenarios_with_main_dimensions)

df_join_scenarios_eigenvalues = scenarios_with_main_dimensions %>% 
  select_at(c("id", scenario_identifier, "sd"))   %>% 
  left_join(
    df_eigenvalues_full,
    by = c("id" = "scenario_id")
  )

# Build data frames that contain eigenvalues and sd
i = 1
while (i <= total_scenarios) {
  current_scenario = scenarios_with_main_dimensions[i, ]
  scenario_id = current_scenario$id
  eigen_data = df_join_scenarios_eigenvalues[df_join_scenarios_eigenvalues$id == scenario_id, ]
  n_main_dimension_scenario = current_scenario$n_main_dimensions
  
  sd_sim = map2(.x = eigen_data$sd, .y = eigen_data$n_main_dimensions, .f = ~ .x[1:.y]) %>% unlist() %>% 
    matrix(., nrow = length(eigen_data$sd), byrow = TRUE) %>% as.data.frame() %>% 
    set_names(nm = paste0("dim_", 1:n_main_dimension_scenario))
  
  sd_sim$scenario_id = scenario_id
  sd_sim$num_sim = eigen_data$num_sim
  sd_sim$method_name = eigen_data$method_name
  
  eigenvalues_sim = map2(.x = eigen_data$eigenvalue_vector, .y = eigen_data$n_main_dimensions, .f = ~ .x[1:.y]) %>% 
    unlist() %>% matrix(., nrow = length(eigen_data$sd), byrow = TRUE) %>% sqrt() %>% as.data.frame() %>%
    set_names(nm = paste0("dim_", 1:n_main_dimension_scenario))
  
  eigenvalues_sim$scenario_id = scenario_id
  eigenvalues_sim$num_sim = eigen_data$num_sim
  eigenvalues_sim$method_name = eigen_data$method_name
  
  if (i == 1) {
    df_sd_flat = sd_sim
    df_eigenvalues_flat = eigenvalues_sim
  } else {
    df_sd_flat = plyr::rbind.fill(df_sd_flat, sd_sim)
    df_eigenvalues_flat = plyr::rbind.fill(df_eigenvalues_flat, eigenvalues_sim)
  }
  
  i = i + 1
}

df_sd_flat_long = df_sd_flat %>% 
  pivot_longer(cols = starts_with("dim"), names_to = "dim", values_to = "sd") %>% 
  filter(!is.na(sd))
 
df_eigenvalues_long = df_eigenvalues_flat %>% 
  pivot_longer(cols = starts_with("dim"), names_to = "dim", values_to = "eigenvalue") %>% 
  filter(!is.na(eigenvalue))


# Enrich scenarios with eigenvalue information
eigenvalues_information = df_join_scenarios_eigenvalues %>% 
  select_at(c("id", scenario_identifier, "num_sim", "method_name")) %>% 
  left_join(df_sd_flat_long, by = c("id" = "scenario_id", "num_sim" = "num_sim", "method_name" = "method_name")) %>% 
  left_join(df_eigenvalues_long, by = c("id" = "scenario_id", "num_sim" = "num_sim", "method_name" = "method_name", "dim" = "dim")) %>% 
  mutate(
    error = eigenvalue - sd
  )

# Analyse data  ----
# MSE
mse_data = eigenvalues_information %>% 
  group_by_at(c(scenario_identifier, "id", "method_name", "dim")) %>% 
  summarise(
    bias_estimator = mean(eigenvalue) - max(sd),
    var_estimator = var(eigenvalue),
    mse_estimator = bias_estimator^2 + var_estimator
  ) %>% 
  arrange_at(c(scenario_identifier, "id"))
  
mse_data %>% View
mse_data %>% ggplot(aes(x = id, color = method_name, y = mse_estimator, group = method_name)) + geom_line()

# Boxplot for the error
i = 1
while (i<=total_scenarios) {
  current_scenario = scenarios_with_main_dimensions[i, ]
  scenario_id = current_scenario$id
  n_main_dimensions = current_scenario$n_main_dimensions
  eigenvalues_information_i = eigenvalues_information %>% filter(id == scenario_id)
  
  name_title = paste0("sample_size: ", current_scenario$sample_size, "; n_cols: ", current_scenario$n_cols, 
                      "; n_main_dim: ", n_main_dimensions, "; sd_main: ", current_scenario$sd_main)
  
  p = eigenvalues_information_i %>% 
    ggplot(aes(x = method_name, y = error, fill = method_name)) + geom_boxplot() + facet_wrap(. ~ dim, ncol = n_main_dimensions) +
    ggtitle(name_title) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none")
  
  png(file.path("plots", "eigenvalues", paste0(name_title, ".png")))
  print(p)
  dev.off()
  i = i + 1
}
