library(tidyverse)
library(ggplot2)
library(ggridges)

# Load data  ----
data_path = file.path(getwd(), 'data')
simulation_data = file.path(data_path, "experiments_processed")
load(file.path(simulation_data, "df_scenarios_full.RData"))
load(file.path(simulation_data, "df_time_full.RData"))
load(file.path(simulation_data, "df_conversion.RData"))
load(file.path(simulation_data, "df_mds_paramenters_full.RData"))

# Manipulate data ----
scenario_identifier = c("sample_size", "n_cols", "n_main_dimensions", "var_main")

# Avoid using scenarions which sample size is 10^6
df_scenarios_full_filtered = df_scenarios_full %>% 
  left_join(df_mds_paramenters_full, by = c("id" = "scenario_id")) %>% 
  filter(!is.na(processed_at), !experiment_label %in% c("l_experiment"))

# Join scenarios and time
df_join_scenarios_time = df_scenarios_full_filtered %>% 
  select_at(c("id", scenario_identifier, "experiment_label", "l_divide", "l_gower", "l_fast")) %>% 
  left_join(
    df_time_full,
    by = c("id" = "scenario_id")
  ) %>% 
  mutate(n_main_dimensions=as.factor(n_main_dimensions),
         n_cols = as.factor(n_cols),
         log_elapsed_time = log10(elapsed_time),
         log_sample_size = log10(sample_size)) %>% 
  left_join(df_conversion)

# Validations
# Each combination of sample_size, n_cols, n_main_dimensions, method_name must appear 100 times
df_join_scenarios_time %>% 
  group_by(sample_size, n_cols, n_main_dimensions, method_name) %>% 
  summarise(total = n()) %>% 
  ungroup() %>% 
  summarise(max_total = max(total), min_total = min(total))

# Plots ----
# Time for l
df_summary_time_l = df_join_scenarios_time %>%
  mutate(x = l_gower) %>% 
  group_by(Algorithm, x) %>%
  summarise(
    mean_elapsed_time = mean(elapsed_time), 
    mean_log_elapsed_time = mean(log_elapsed_time)
  )
  #) %>% mutate(x = as.factor(x))

levels(df_summary_time_l$Algorithm) <- c("D&C", "Interp", "Fast", "RMDS")

df_summary_time_l %>% 
  ggplot(aes(x = x, y = mean_elapsed_time, group = Algorithm, color = Algorithm)) +
  geom_point(size = 2) +
  geom_line() +
  theme(panel.spacing.y=unit(0.5, "lines"), legend.position="bottom") +
  xlab("\u2113 value") + 
  ylab("Elapsed time (sec.)") +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91")) +
  xlim(c(0, 1600))
  #+ ggsave(file.path(getwd(), "images", "l_parameter_time.png"), dev = 'png', width = 10, height = 10, units = "cm")

# Time for sample size
df_summary_time_sample_size = df_join_scenarios_time %>%
  mutate(x = log_sample_size) %>% 
  group_by(Algorithm, x) %>%
  summarise(
    mean_elapsed_time = mean(elapsed_time), 
    mean_log_elapsed_time = mean(log_elapsed_time)
  ) 

levels(df_summary_time_sample_size$Algorithm) <- c("D&C", "Interpolation MDS", "Fast", "RMDS (Paradis 2021)")

df_summary_time_sample_size %>% 
  ggplot(aes(x = x, y = mean_log_elapsed_time, group = Algorithm, color = Algorithm)) +
  geom_point(size = 2) +
  geom_line() +
  theme(panel.spacing.y=unit(0.5, "lines"), legend.position="bottom") +
  xlab("log10 of sample size") + 
  ylab("Mean of log10 of elapsed time (seconds)") +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91"))
  #ggsave(file.path(getwd(), "images",  "mean_log_time_all.png"), dev = 'png', width = 22, height = 18, units = "cm")

#Log scale
df_summary_time_sample_size_log_scale = df_join_scenarios_time %>%
  mutate(x = sample_size) %>% 
  group_by(Algorithm, x) %>%
  summarise(
    mean_elapsed_time = mean(elapsed_time), 
    mean_log_elapsed_time = mean(log_elapsed_time)
  ) 

levels(df_summary_time_sample_size_log_scale$Algorithm) <- c("D&C", "Interpolation MDS", "Fast", "RMDS (Paradis 2021)")


df_summary_time_sample_size_log_scale %>% 
  ggplot(aes(x = x, y = mean_elapsed_time, group = Algorithm, color = Algorithm)) +
  geom_point(size = 2) +
  geom_line() +
  theme(panel.spacing.y=unit(0.5, "lines"), legend.position="bottom") +
  scale_x_log10() +
  scale_y_log10() + 
  xlab("sample size") + 
  ylab("Mean of elapsed time (in seconds)") +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91"))
ggsave(file.path(getwd(), "images",  "mean_time_all_log_scale.png"), dev = 'png', width = 22, height = 18, units = "cm")


# Quantiles: Particular case
quant = c(0.025, 1-0.025)
quant_char = c(paste0("q_", as.character(quant)), "mean")
df_join_scenarios_time %>% 
  filter(sample_size == 10^6, n_cols == 100, n_main_dimensions == 10) %>% 
  group_by(Algorithm) %>% 
  summarise(
    statistic = quant_char,
    value = c(quantile(elapsed_time, quant), mean(elapsed_time))) %>% 
  pivot_wider(
    names_from = statistic,
    values_from = value)
