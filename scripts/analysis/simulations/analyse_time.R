library(tidyverse)
library(ggplot2)
library(ggridges)

# Load data  ----
data_path = file.path(getwd(), 'data')
load(file.path(data_path, "df_scenarios_full.RData"))
load(file.path(data_path, "df_time_full.RData"))
load(file.path(data_path, "df_conversion.RData"))
load(file.path(data_path, "df_mds_paramenters_full.RData"))

# Manipulate data ----
scenario_identifier = c("sample_size", "n_cols", "n_main_dimensions", "var_main")

# Avoid using scenarions which sample size is 10^6
df_scenarios_full_filtered = df_scenarios_full %>% 
  left_join(df_mds_paramenters_full, by = c("id" = "scenario_id")) %>% 
  filter(!is.na(processed_at), experiment_label %in% c("l_experiment"))

# Join scenarios and time
df_join_scenarios_time = df_scenarios_full_filtered %>% 
  select_at(c("id", scenario_identifier, "experiment_label", "l")) %>% 
  left_join(
    df_time_full,
    by = c("id" = "scenario_id")
  ) %>% 
  mutate(n_main_dimensions=as.factor(n_main_dimensions),
         n_cols = as.factor(n_cols),
         log_elapsed_time = log(elapsed_time)) %>% 
  left_join(df_conversion)



# Analyse data  ----
# Plots by sample size
df_summary_sample_size = df_join_scenarios_time %>%
  group_by(sample_size, algorithm) %>%
  summarise(
    mean_elapsed_time = mean(elapsed_time), 
    mean_log_elapsed_time = mean(log_elapsed_time)
  )


df_summary_sample_size %>% View()

df_summary_sample_size %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  ggplot(aes(x = sample_size, y = mean_log_elapsed_time, group = algorithm, color = algorithm)) +
  geom_point(size = 2) +
  geom_line() +
  theme(panel.spacing.y=unit(0.5, "lines"), legend.position="bottom") +
  xlab("Sample Size") + 
  ylab("Mean of Log. Elapsed time") +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91")) +
  ggsave(file.path(getwd(), "images", "mean_log_time_all.png"),
         dpi = 300, dev = 'png', height = 18, width = 22, units = "cm")

# Quantiles
quant = c(0.025, 1-0.025)
quant_char = c(paste0("q_", as.character(quant)), "mean")
df_join_scenarios_time %>% 
  filter(sample_size == 10^6, n_cols == 100, n_main_dimensions == 10) %>% 
  group_by(algorithm) %>% 
  summarise(
    statistic = quant_char,
    value = c(quantile(elapsed_time, quant), mean(elapsed_time))) %>% 
  pivot_wider(
    names_from = statistic,
    values_from = value)

df_join_scenarios_time %>% 
  filter(sample_size == 10^5, n_cols == 100, n_main_dimensions == 10) %>% 
  group_by(algorithm) %>% 
  summarise(
    statistic = quant_char,
    value = c(quantile(elapsed_time, quant), mean(elapsed_time))) %>% 
  pivot_wider(
    names_from = statistic,
    values_from = value)
