library(tidyverse)
library(ggplot2)

# Load data ----
data_path = file.path(getwd(), 'data')
simulation_data = file.path(data_path, "experiments_processed")
load(file.path(simulation_data, "df_scenarios_full.RData"))
load(file.path(simulation_data, "df_eigenvalues_full.RData"))
load(file.path(simulation_data, "df_conversion.RData"))
load(file.path(simulation_data, "df_mds_paramenters_full.RData"))

# Manipulate data ----
scenario_identifier = c("sample_size", "n_cols", "n_main_dimensions", "var_main")

# Join scenarios and eigenvalues
scenarios_with_main_dimensions = df_scenarios_full %>% 
  left_join(df_mds_paramenters_full, by = c("id" = "scenario_id")) %>% 
  filter(n_main_dimensions > 0,  !is.na(processed_at), !experiment_label %in% c("l_experiment"))

total_scenarios = nrow(scenarios_with_main_dimensions)

df_join_scenarios_eigenvalues = scenarios_with_main_dimensions %>% 
  select_at(c("id", scenario_identifier, "var", "experiment_label"))   %>% 
  left_join(
    df_eigenvalues_full,
    by = c("id" = "scenario_id")) 

# Build data frames that contain eigenvalues and var
i = 1
while (i <= total_scenarios) {
  current_scenario = scenarios_with_main_dimensions[i, ]
  scenario_id = current_scenario$id
  eigen_data = df_join_scenarios_eigenvalues[df_join_scenarios_eigenvalues$id == scenario_id, ]
  n_main_dimension_scenario = current_scenario$n_main_dimensions
  
  var_sim = map2(.x = eigen_data$var, .y = eigen_data$n_main_dimensions, .f = ~ .x[1:.y]) %>% unlist() %>% 
    matrix(., nrow = length(eigen_data$var), byrow = TRUE) %>% as.data.frame() %>% 
    set_names(nm = paste0("dim_", sprintf("%02d", 1:n_main_dimension_scenario)))
  
  var_sim$scenario_id = scenario_id
  var_sim$num_sim = eigen_data$num_sim
  var_sim$method_name = eigen_data$method_name
  
  eigenvalues_sim = map2(.x = eigen_data$eigenvalue_vector, .y = eigen_data$n_main_dimensions, .f = ~ .x[1:.y]) %>% 
    unlist() %>% matrix(., nrow = length(eigen_data$var), byrow = TRUE) %>% as.data.frame() %>%
    set_names(nm = paste0("dim_", sprintf("%02d", 1:n_main_dimension_scenario)))
  
  eigenvalues_sim$scenario_id = scenario_id
  eigenvalues_sim$num_sim = eigen_data$num_sim
  eigenvalues_sim$method_name = eigen_data$method_name
  
  if (i == 1) {
    df_var_flat = var_sim
    df_eigenvalues_flat = eigenvalues_sim
  } else {
    df_var_flat = plyr::rbind.fill(df_var_flat, var_sim)
    df_eigenvalues_flat = plyr::rbind.fill(df_eigenvalues_flat, eigenvalues_sim)
  }
  
  i = i + 1
}

df_var_flat_long = df_var_flat %>% 
  pivot_longer(cols = starts_with("dim"), names_to = "dim", values_to = "var") %>% 
  filter(!is.na(var))
 
df_eigenvalues_long = df_eigenvalues_flat %>% 
  pivot_longer(cols = starts_with("dim"), names_to = "dim", values_to = "eigenvalue") %>% 
  filter(!is.na(eigenvalue))

# Enrich scenarios with eigenvalue information
eigenvalues_information = df_join_scenarios_eigenvalues %>% 
  select_at(c("id", scenario_identifier, "num_sim", "method_name", "experiment_label")) %>% 
  left_join(df_var_flat_long, by = c("id" = "scenario_id", "num_sim" = "num_sim", "method_name" = "method_name")) %>% 
  left_join(df_eigenvalues_long, by = c("id" = "scenario_id", "num_sim" = "num_sim", "method_name" = "method_name", "dim" = "dim")) %>% 
  left_join(df_conversion) %>% 
  left_join(df_mds_paramenters_full, by = c("id" = "scenario_id")) %>% 
  mutate(error = eigenvalue - var,
         error_2 = error^2)

levels(eigenvalues_information$Algorithm) <- c("D&C", "Interp", "Fast", "RMDS")


# Plots  ----
# RMSE Eigenvalues for l parameter
eigenvalues_information %>%
  group_by(Algorithm, l_divide)%>%
  summarise(rmse = sqrt(mean(error_2))) %>% 
  mutate(x = l_divide) %>% 
  #mutate(x = as.factor(x)) %>% 
  ggplot(aes(x = x, y = rmse, group = Algorithm, color = Algorithm)) +
  geom_line() +
  geom_point(size = 2) +
  theme(panel.spacing.y=unit(0.5, "lines"), legend.position="bottom") +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91")) + 
  xlab("\u2113 value") + 
  ylab("RMSE Eigenvalues") +
  xlim(c(0, 1600))
#+ggsave(file.path(getwd(), "images", "l_parameter_bias.png"), dev = 'png', width = 10, height = 10, units = "cm")
 
# Bias
eigenvalues_information %>% 
  group_by(sample_size, Algorithm, dim) %>% 
  summarise(bias = mean(error)) %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  ggplot(aes(x = sample_size, y = bias, group = Algorithm, color = Algorithm)) +
  geom_point() +
  facet_wrap( ~ dim, ncol = 2) +
  theme(panel.spacing.y=unit(0.5, "lines"), legend.position = "bottom") +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91")) + 
  xlab("Sample Size") + 
  ylab("Bias") +
  geom_hline(yintercept = 0, linetype = "dashed") 
  #ggsave(file.path(getwd(), "images", "bias_all.png"),dpi = 300, dev = 'png', height = 18, width = 22, units = "cm")

# Using another variable in the facet_wrap
#pdf('images/bias_all_by_n.pdf', width = 4, height = 8)
eigenvalues_information %>% 
  group_by(sample_size, Algorithm, dim) %>% 
  summarise(bias = mean(error)) %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  ggplot(aes(x = dim, y = bias, group = Algorithm, color = Algorithm)) +
  geom_point() +
  facet_wrap( ~ sample_size, ncol = 1) +
  theme(
    panel.spacing.y=unit(0.5, "lines"), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91")) + 
  xlab("Dimension") + 
  ylab("Bias") +
  geom_hline(yintercept = 0, linetype = "dashed")
#dev.off()

# MSE/RMSE
eigenvalues_information %>% 
  group_by(sample_size, Algorithm, dim) %>% 
  summarise(
    mse = mean(error^2),
    rmse = sqrt(mse)
  ) %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  ggplot(aes(x = sample_size, y = rmse, group = Algorithm, color = Algorithm)) +
  geom_point() +
  facet_wrap( ~ dim, ncol = 2) +
  theme(panel.spacing.y=unit(0.5, "lines"), legend.position="bottom") +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91")) + 
  xlab("Dimension") + 
  ylab("RMSE")
  #ggsave(file.path(getwd(), "images", "rmse_all.png"),dpi = 300, dev = 'png', height = 18, width = 22, units = "cm")


# Using another variable in the facet_wrap
#pdf('images/rmse_all_by_n.pdf', width = 4, height = 8)
eigenvalues_information %>% 
  group_by(sample_size, Algorithm, dim) %>% 
  summarise(
    mse = mean(error^2),
    rmse = sqrt(mse)
  ) %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  ggplot(aes(x = dim, y = rmse, group = Algorithm, color = Algorithm)) +
  geom_point() +
  facet_wrap( ~ sample_size, ncol = 1) +
  theme(
    panel.spacing.y=unit(0.5, "lines"),
    legend.position="bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91")) + 
  xlab("Dimension") + 
  ylab("RMSE")
#dev.off()

# Error for a particular scenario
eigenvalues_information %>% 
  filter(sample_size == 10^6, n_main_dimensions == 10, n_cols == 100) %>% 
  ggplot(aes(x = Algorithm, y = error, fill = Algorithm)) +
  geom_boxplot() + 
  facet_wrap( ~ dim, ncol = 5) +
  theme(
    panel.spacing.y=unit(0.5, "lines"), 
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.position="bottom"
  ) +
  scale_fill_manual(values = c("#0000FF", "#FF0000", "#00AF91")) + 
  ylab("Error") +
  geom_hline(yintercept = 0, linetype = "dashed")
  #ggsave(file.path(getwd(), "images", "error_1000000_100_10.png"),dpi = 300, dev = 'png', height = 18, width = 22, units = "cm")
