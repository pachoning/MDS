options(digits = 12)   
library(tidyverse)
library(ggplot2)

# Load data ----
data_path = file.path(getwd(), 'data')
load(file.path(data_path, "df_scenarios_full.RData"))
load(file.path(data_path, "df_correlation_full.RData"))
load(file.path(data_path, "df_conversion.RData"))
load(file.path(data_path, "df_mds_paramenters_full.RData"))

# Manipulate data ----
scenario_identifier = c("sample_size", "n_cols", "n_main_dimensions", "var_main")

# Avoid using scenarions which sample size is 10^6
scenarios_with_main_dimesions = df_scenarios_full %>% 
  left_join(df_mds_paramenters_full, by = c("id" = "scenario_id")) %>% 
  filter(n_main_dimensions > 0, !is.na(processed_at), experiment_label %in% c("l_experiment"))

# Join scenarios and correlation
df_join_scenarios_correlation = scenarios_with_main_dimesions %>% 
  select_at(c("id", scenario_identifier, "experiment_label", "l")) %>% 
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
  group_by(sample_size, n_cols, n_main_dimensions, algorithm, dim) %>% 
  summarise(total = n()) %>% 
  filter(total > 100)


# Analyse data  ----
# Table of correlations
quant = c(0.025, 1-0.025)
quant_char = c(paste0("q_", as.character(quant)), "mean")
df_join_scenarios_correlation %>% 
  group_by(algorithm) %>% 
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

# Plots
# Plot of the correlation
df_join_scenarios_correlation %>% 
  group_by(algorithm, sample_size, dim) %>%
  summarise(mean_correlation = mean(correlation)) %>%
  mutate(sample_size = as.factor(sample_size)) %>% 
  ggplot(aes(x = sample_size, y = mean_correlation, color = algorithm, group = algorithm)) +
  geom_point() +
  facet_wrap( ~ dim, ncol = 2) +
  theme(panel.spacing.y = unit(0.5, "lines"), legend.position="bottom") +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91")) + 
  xlab("Sample Size") + 
  ylab("Mean Correlation") +
  ggsave(file.path(getwd(), "images", "correlation_all.png"),
         dpi = 300, dev = 'png', height = 18, width = 22, units = "cm")
  
# Particular cases
df_join_scenarios_correlation %>% 
  filter(sample_size == 10^6, n_main_dimensions == 10, n_cols == 100) %>% 
  ggplot(aes(x = algorithm, y = correlation, fill = algorithm)) + 
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
