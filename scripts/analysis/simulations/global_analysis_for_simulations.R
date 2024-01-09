library(tidyverse)
library(xtable)

paths_to_data <- c(
  #file.path(getwd(), "data", "full_experiments")
  #file.path(getwd(), "data", "full_experiments_new_l_10"),
  #file.path(getwd(), "data", "full_experiments_new_l_20"),
  #file.path(getwd(), "data", "full_experiments_new_l_30"),
  #file.path(getwd(), "data", "full_experiments_new_l_40")
  #file.path(getwd(), "data", "full_experiments_80")
  #file.path(getwd(), "data", "full_experiments_15000"),
  #file.path(getwd(), "data", "full_experiments_25000_to_95000")
  file.path(getwd(), "data", "full_experiments_v2_10"),
  file.path(getwd(), "data", "full_experiments_v2_90"),
  file.path(getwd(), "data", "full_experiments_v2_250_500_750")
)

order_algorithms <- c(
  "landmark_mds",
  "interpolation_mds",
  "reduced_mds",
  "pivot_mds",
  "divide_conquer_mds",
  "fast_mds"
)

i_paths <- 1
for (path in paths_to_data) {
  load(file.path(path, "correlations.RData"))
  load(file.path(path, "eigenvalues.RData"))
  load(file.path(path, "scenarios.RData"))
  load(file.path(path, "times.RData"))
  load(file.path(path, "params.RData"))
  load(file.path(path, "l_param.RData"))
  if(i_paths == 1) {
    df_correlations <- correlations
    df_eigenvalues <- eigenvalues
    df_scenarios <- scenarios
    df_times <- times
    df_params <- params
    df_l_param <- l_param
    
  } else{
    df_correlations <- rbind(df_correlations, correlations)
    df_eigenvalues <- rbind(df_eigenvalues, eigenvalues)
    df_scenarios <- rbind(df_scenarios, scenarios)
    df_times <- rbind(df_times, times)
    df_params <- rbind(df_params, params)
    df_l_param <- rbind(df_l_param, l_param)
  }
  print(nrow(df_times))
  i_paths <- i_paths + 1
}

# Sanity check
if (length(unique(df_scenarios$id)) != nrow(df_scenarios)) {
  vector_repeated_value <- df_scenarios$id[duplicated(df_scenarios$id)]
  repeated_value <- paste0(vector_repeated_value, collapse = ", ")
  stop(paste0("Repeated values for scenario_id: ", repeated_value))
}

# Visualise l_param
df_l_param %>% 
  group_by(algorithm) %>% 
  summarise(
    min_val = min(l_value),
    mean_val = mean(l_value),
    max_val = max(l_value)
  ) %>% 
  View()

df_scenarios %>% 
  distinct(n_sample, variance) %>% 
  View()

# Analysis of correlation ----
df_scenario_corr <- df_scenarios %>% 
  left_join(df_correlations, by = c("id" = "scenario_id")) %>% 
  mutate(
    algorithm = factor(algorithm, levels = c(order_algorithms))
  ) %>% 
  unnest(correlation_coeffs) %>% 
  group_by(id, algorithm, simulation_id) %>% 
  mutate(
    n_sample = unlist(n_sample),
    n_dim = lengths(variance),
    main_dim = 1:n(),
    n_main_dim = map_int(variance, function(x) sum(x != 1))
  ) %>% 
  ungroup()

View(df_scenario_corr)

# Sanity checks
df_sanity_check <- df_scenario_corr %>% 
  distinct(id, simulation_id, algorithm, n_sample, n_main_dim, n_dim) %>% 
  group_by(algorithm, n_sample, n_main_dim, n_dim) %>% 
  summarise(
    total = n(),
    min_total = min(total),
    max_total = max(total)
  ) %>% 
  ungroup() 

df_sanity_check %>% View()

df_sanity_check %>% 
  distinct(n_sample, n_main_dim, n_dim) %>% 
  View()

df_sanity_check %>% 
  distinct(min_total, max_total)

# Correlation by algorithm
df_scenario_corr %>% 
  group_by(algorithm) %>% 
  summarise(
    q1_corr = quantile(correlation_coeffs, 0.025),
    mean_corr = mean(correlation_coeffs),
    q3 = quantile(correlation_coeffs, 1 - 0.025)
  ) %>% 
  xtable(digits = 5)


# Count the number of correlation coefficients for a single interation of all scenarios
df_scenario_corr %>% 
  distinct(algorithm, n_sample, n_dim, n_main_dim, main_dim) %>%
  group_by(algorithm) %>% 
  count()

# Count the number of correlation coefficients taking into account all scenarios and all simulations
df_scenario_corr %>% 
  distinct(algorithm, n_sample, n_dim, n_main_dim, main_dim, simulation_id, id) %>%
  group_by(algorithm) %>% 
  count() 

# Table for small sample sizes
sample_size_filter <- 100000
df_scenario_corr %>% 
  #filter(n_sample <= sample_size_filter) %>% 
  group_by(n_sample, algorithm) %>% 
  summarise(
    mean_corr_coef = mean(correlation_coeffs),
    min_corr_coef = min(correlation_coeffs)
  ) %>% View

# Plot for small sample sizes
df_scenario_corr %>% 
  filter(n_sample <= sample_size_filter) %>% 
  ggplot(aes(x = n_sample, y = -log(1-correlation_coeffs), group = interaction(n_sample, algorithm))) + 
  geom_boxplot(
    aes(colour = algorithm),
    width = 3000
  ) + 
  stat_summary(
    aes(group = algorithm, colour = algorithm), 
    fun = median,
    geom="line",
    position = position_dodge(width = 3000)
  )

# Plot for small sample sizes splitted by main dimensions
df_scenario_corr %>% 
  filter(n_sample <= sample_size_filter) %>% 
  ggplot(aes(x = n_sample, y = -log(1-correlation_coeffs), group = interaction(n_sample, algorithm))) + 
  geom_boxplot(
    aes(colour = algorithm),
    width = 3000
  ) + 
  stat_summary(
    aes(group = algorithm, colour = algorithm), 
    fun = median,
    geom="line",
    position = position_dodge(width = 3000)
  ) +
  facet_wrap(vars(n_dim, n_main_dim), labeller = "label_both") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Plot for large sample sizes splitted by main dimensions
df_scenario_corr %>% 
  filter(n_sample > sample_size_filter) %>% 
  ggplot(aes(x = n_sample, y = -log(1-correlation_coeffs), group = interaction(n_sample, algorithm))) + 
  geom_boxplot(
    aes(colour = algorithm),
    width = 3000
  ) + 
  stat_summary(
    aes(group = algorithm, colour = algorithm), 
    fun = median,
    geom="line",
    position = position_dodge(width = 3000)
  ) +
  facet_wrap(vars(n_dim, n_main_dim), labeller = "label_both") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# Analysis of time ----
df_scenario_times <- df_scenarios %>% 
  left_join(df_times, by = c("id" = "scenario_id")) %>% 
  mutate(
    algorithm = factor(algorithm, levels = c(order_algorithms)),
    n_sample = unlist(n_sample),
    n_dim = lengths(variance),
    n_main_dim = map_int(variance, function(x) sum(x != 1))
  )
View(df_scenario_times)

# Sanity checks
df_scenario_times %>% 
  distinct(id, simulation_id, algorithm, n_sample, n_main_dim, n_dim) %>% 
  group_by(algorithm, n_sample, n_main_dim, n_dim) %>% 
  summarise(
    total = n(),
    min_total = min(total),
    max_total = max(total)
  ) %>% 
  ungroup() %>% 
  distinct(min_total, max_total)

# Plot for time by dominant dimensions
pdf('images/mean_time_all_log_scale_6.pdf', width = 6, height = 6)
df_scenario_times %>% 
  group_by(algorithm, n_sample) %>%
  summarise(
    mean_elapsed_time = mean(elapsed_time), 
    mean_log_elapsed_time = mean(elapsed_time)
  ) %>% 
  ggplot(aes(x = n_sample, y = mean_elapsed_time, group = algorithm, color = algorithm)) +
  geom_point() +
  geom_line() +
  theme(panel.spacing.y=unit(0.5, "lines"), legend.position="bottom") +
  scale_x_log10() +
  scale_y_log10() + 
  xlab("sample size") + 
  ylab("Mean of elapsed time (in seconds)")
dev.off()

df_scenario_times %>% 
  group_by(algorithm) %>% 
  summarise(
    value = mean(elapsed_time/n_sample) * 10^6
  ) %>% View()

df_scenario_times %>% 
  filter(n_sample == 1000000, n_dim == 100, n_main_dim == 10) %>% 
  group_by(algorithm) %>% 
  summarise(
    q1 = quantile(elapsed_time, 0.025),
    mean_val = mean(elapsed_time),
    q3 = quantile(elapsed_time, 1 - 0.025)
  ) %>% 
    xtable(digits = 3)

df_scenario_times %>% 
  mutate(var_sce = paste0(sprintf("%03d", n_dim), "--", sprintf("%02d", n_main_dim))) %>% 
  filter(n_sample == 1000000) %>% 
  arrange(var_sce) %>% 
  group_by(algorithm) %>% 
  mutate(mean_all = mean(elapsed_time)) %>% 
  ungroup() %>% 
  group_by(algorithm, var_sce) %>% 
  summarise(
    mean_val = mean(elapsed_time),
    mean_all = max(mean_all)
  ) %>%
  pivot_wider(names_from = var_sce, values_from = c(mean_val, mean_all)) %>% View()
  

# Plot for small sample sizes splitted by main dimensions
df_scenario_times %>% 
  filter(n_sample <= sample_size_filter) %>% 
  ggplot(aes(x = n_sample, y = elapsed_time, group = interaction(n_sample, algorithm))) + 
  geom_boxplot(
    aes(colour = algorithm),
    width = 3000
  ) + 
  stat_summary(
    aes(group = algorithm, colour = algorithm), 
    fun = median,
    geom="line",
    position = position_dodge(width = 3000)
  ) +
  facet_wrap(vars(n_dim, n_main_dim), labeller = "label_both") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Plot for small sample sizes splitted by main dimensions without fast and pivot
df_scenario_times %>% 
  filter(n_sample <= sample_size_filter, !algorithm %in% c("fast_mds", "pivot_mds")) %>% 
  ggplot(aes(x = n_sample, y = elapsed_time, group = interaction(n_sample, algorithm))) + 
  geom_boxplot(
    aes(colour = algorithm),
    width = 3000
  ) + 
  stat_summary(
    aes(group = algorithm, colour = algorithm), 
    fun = median,
    geom="line",
    position = position_dodge(width = 3000)
  ) +
  facet_wrap(vars(n_dim, n_main_dim), labeller = "label_both") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

df_scenario_times %>% 
  filter(n_sample > sample_size_filter) %>% 
  ggplot(aes(x = n_sample, y = elapsed_time, group = interaction(n_sample, algorithm))) + 
  geom_boxplot(
    aes(colour = algorithm),
    width = 3000
  ) + 
  stat_summary(
    aes(group = algorithm, colour = algorithm), 
    fun = median,
    geom="line",
    position = position_dodge(width = 3000)
  ) +
  facet_wrap(vars(n_dim, n_main_dim), labeller = "label_both") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Analysis of eigenvalues ----
# According to Wikipedia (and LyD), bias = E(estimator) - parameter
df_scenario_eigen <- df_scenarios %>% 
  left_join(df_eigenvalues, by = c("id" = "scenario_id")) %>% 
  mutate(
    algorithm = factor(algorithm, levels = c(order_algorithms))
  ) %>% 
  unnest(eigenvalues) %>% 
  group_by(id, algorithm, simulation_id) %>% 
  mutate(
    n_sample = unlist(n_sample),
    n_dim = lengths(variance),
    main_dim = 1:n(),
    n_main_dim = map_int(variance, function(x) sum(x != 1))
  ) %>% 
  ungroup()

View(df_scenario_eigen)

# Sanity checks
df_scenario_eigen %>% 
  distinct(id, simulation_id, algorithm, n_sample, n_main_dim, n_dim) %>% 
  group_by(algorithm, n_sample, n_main_dim, n_dim) %>% 
  summarise(
    total = n(),
    min_total = min(total),
    max_total = max(total)
  ) %>% 
  ungroup() %>% 
  distinct(min_total, max_total)

# Plot bias by dimensions
#pdf(
#  "/Users/cristianpachon/MEGA/tesis_cristian/mds_for_big_data/to_ADAC/bias_all_by_n_6_line.pdf",
#  width = 5,
#  height = 8
#)
pdf('images/bias_all_by_n_6_line.pdf', width = 5, height = 8)
df_scenario_eigen %>% 
  group_by(n_sample, algorithm, main_dim, n_main_dim) %>% 
  summarise(bias = mean(eigenvalues) - 15) %>% 
  mutate(
    sample_size = as.factor(n_sample),
    main_dim = as.factor(main_dim)
  ) %>% 
  ggplot(aes(x = main_dim, y = bias, group = algorithm, color = algorithm)) +
  geom_point() +
  geom_line() + 
  facet_grid( 
    cols = vars(n_main_dim),
    rows = vars(sample_size),
    scales = "free_x", 
    space = 'free'
  ) +
  theme(
    panel.spacing.y=unit(0.5, "lines"), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  xlab("Dimension") + 
  ylab("Bias") +
  geom_hline(yintercept = 0, linetype = "dashed") 
dev.off()

# Plot RMSE by dimensions
pdf('images/rmse_all_by_n_6_line.pdf', width = 5, height = 8)
df_scenario_eigen %>% 
  group_by(n_sample, algorithm, main_dim, n_main_dim) %>% 
  summarise(
    bias = mean(eigenvalues) - 15,
    rmse = sqrt(var(eigenvalues) + bias^2)
  ) %>% 
  mutate(
    sample_size = as.factor(n_sample),
    main_dim = as.factor(main_dim)
  ) %>% 
  ggplot(aes(x = main_dim, y = rmse, group = algorithm, color = algorithm)) +
  geom_point() +
  geom_line() + 
  facet_grid( 
    cols = vars(n_main_dim),
    rows = vars(sample_size),
    scales = "free_x", 
    space = 'free'
  ) +
  theme(
    panel.spacing.y=unit(0.5, "lines"), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  xlab("Dimension") + 
  ylab("RMSE") +
  geom_hline(yintercept = 0, linetype = "dashed") 
dev.off()


# Plot for small sample sizes splitted by main dimensions
df_scenario_eigen %>% 
  filter(n_sample <= sample_size_filter) %>% 
  ggplot(aes(x = n_sample, y = bias, group = interaction(n_sample, algorithm))) + 
  geom_boxplot(
    aes(colour = algorithm),
    width = 3000
  ) + 
  stat_summary(
    aes(group = algorithm, colour = algorithm), 
    fun = median,
    geom="line",
    position = position_dodge(width = 3000)
  ) +
  facet_wrap(vars(n_dim, n_main_dim), labeller = "label_both") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
