options(digits = 12)   
library(tidyverse)
library(ggplot2)

# Load data ----
data_path = file.path(getwd(), 'data')
simulation_data = file.path(data_path, "experiments_processed")
load(file.path(simulation_data, "df_scenarios_full.RData"))
load(file.path(simulation_data, "df_correlation_full.RData"))
load(file.path(simulation_data, "df_conversion.RData"))
load(file.path(simulation_data, "df_mds_paramenters_full.RData"))

# Manipulate data ----
scenario_identifier = c("sample_size", "n_cols", "n_main_dimensions", "var_main")

# Avoid using scenarions which sample size is 10^6
scenarios_with_main_dimesions = df_scenarios_full %>% 
  left_join(df_mds_paramenters_full, by = c("id" = "scenario_id")) %>% 
  filter(n_main_dimensions > 0, !is.na(processed_at), experiment_label %in% c("l_experiment"))

# Join scenarios and correlation
df_join_scenarios_correlation = scenarios_with_main_dimesions %>% 
  select_at(c("id", scenario_identifier, "experiment_label", "l_divide", "l_gower", "l_fast")) %>% 
  left_join(
    df_correlation_full %>% select(-n_main_dimensions),
    by = c("id" = "scenario_id")
  ) %>% unnest(correlation_vector) %>% 
  rename(correlation = correlation_vector) %>% 
  group_by(id, method_name, num_sim) %>% 
  mutate(dim = paste0("dim_", sprintf("%02d", 1:n()))) %>% 
  left_join(df_conversion)

# Validations ----
# Unique combinations
df_join_scenarios_correlation %>% 
  group_by(sample_size, n_cols, n_main_dimensions, Algorithm, dim) %>% 
  summarise(total = n()) %>% 
  ungroup() %>% 
  summarise(min_total = min(total), max_total = max(total))

levels(df_join_scenarios_correlation$Algorithm) <- c("D&C", "Interp", "Fast", "RMDS")

# Analyse data  ----
# Table of correlations
quant = c(0.025, 1-0.025)
quant_char = c(paste0("q_", as.character(quant)), "mean")
df_join_scenarios_correlation %>% 
  group_by(Algorithm) %>% 
  summarise(
    statistic = quant_char,
    value = c(quantile(correlation, quant), mean(correlation))) %>% 
  pivot_wider(
    names_from = statistic,
    values_from = value) %>% 
  mutate_if(
    is.numeric,
    function(x) formatC(x, format = "e", digits = 5)
  )

# Plots ----
# Plot of the correlation for l
df_join_scenarios_correlation %>% 
  group_by(Algorithm, l_divide) %>%
  summarise(mean_correlation = mean(correlation)) %>%
  mutate(x = l_divide) %>%
  #mutate(x = as.factor(x))%>% 
  ggplot(aes(x = x, y = mean_correlation, color = Algorithm, group = Algorithm)) +
  geom_line() +
  geom_point(size = 2) +
  theme(panel.spacing.y = unit(0.5, "lines"), legend.position="bottom") +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91")) + 
  xlab("\u2113 value") + 
  ylab("Correlation") +
  xlim(c(0, 1600)) +
  ylim(c(0.888, 1)) 
#+ggsave(file.path(getwd(), "images", "l_parameter_correlation.png"), dev = 'png', width = 10, height = 10, units = "cm")

# Plot of the correlation
df_join_scenarios_correlation %>% 
  group_by(Algorithm, sample_size, dim) %>%
  summarise(mean_correlation = mean(correlation)) %>%
  mutate(x = sample_size) %>%
  mutate(x = as.factor(x)) %>% 
  ggplot(aes(x = x, y = mean_correlation, color = Algorithm, group = Algorithm)) +
  geom_point() +
  facet_wrap( ~ dim, ncol = 2) +
  theme(panel.spacing.y = unit(0.5, "lines"), legend.position="bottom") +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91")) + 
  xlab("Sample Size") + 
  ylab("Correlation") +
  ggsave(file.path(getwd(), "images", "correlation_all.png"),  dev = 'png', height = 18, width = 22, units = "cm")
  
# Particular cases
df_join_scenarios_correlation %>% 
  filter(sample_size == 10^6, n_main_dimensions == 10, n_cols == 100) %>% 
  ggplot(aes(x = Algorithm, y = correlation, fill = Algorithm)) + 
  geom_boxplot() + 
  facet_wrap(. ~ dim, ncol = 5) +
  theme(
    panel.spacing.y = unit(0.5, "lines"), 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom") +
  ylab("Correlation coefficient") +
  scale_fill_manual(values = c("#0000FF", "#FF0000", "#00AF91")) + 
  ggsave(file.path(getwd(), "images", "correlation_1000000_100_10.png"),
         dpi = 300, dev = 'png', height = 18, width = 22, units = "cm")
