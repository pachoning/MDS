library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)

paths_to_data <- c(
  #file.path(getwd(), "data", "l_param_all"),
  #file.path(getwd(), "data", "l_param_900"),
  #file.path(getwd(), "data", "l_param_1000"),
  #file.path(getwd(), "data", "l_param_1100"),
  #file.path(getwd(), "data", "l_param_50_to_800")
  #file.path(getwd(), "data", "l_param_100_to_400_factor_5"),
  #file.path(getwd(), "data", "l_param_100_to_400_factor_5_v2"
  file.path(getwd(), "data", "l_param_100_to_1000_factor_5"),
  file.path(getwd(), "data", "l_param_100_to_1000_factor_5_19")
  #file.path(getwd(), "data", "l_param_800_others_factor_5_20")
)

order_algorithms <- c(
  "landmark_mds",
  "interpolation_mds",
  "reduced_mds",
  "pivot_mds",
  "divide_conquer_mds",
  "fast_mds"
)

# Join all the information
i_paths <- 1
for (path in paths_to_data) {
  load(file.path(path, "correlations.RData"))
  load(file.path(path, "eigenvalues.RData"))
  load(file.path(path, "scenarios.RData"))
  load(file.path(path, "times.RData"))
  load(file.path(path, "params.RData"))
  
  if(i_paths == 1) {
    df_correlations <- correlations
    df_eigenvalues <- eigenvalues
    df_scenarios <- scenarios
    df_times <- times
  } else{
    df_correlations <- rbind(df_correlations, correlations)
    df_eigenvalues <- rbind(df_eigenvalues, eigenvalues)
    df_scenarios <- rbind(df_scenarios, scenarios)
    df_times <- rbind(df_times, times)
  }
  i_paths <- i_paths + 1
}

# Sanity check
if (length(unique(df_scenarios$id)) != nrow(df_scenarios)) {
  vector_repeated_value <- df_scenarios$id[duplicated(df_scenarios$id)]
  repeated_value <- paste0(vector_repeated_value, collapse = ", ")
  stop(paste0("Repeated values for scenario_id: ", repeated_value))
}

# Analyse correlation ----
df_scenario_corr <- df_scenarios %>% 
  left_join(df_correlations, by = c("id" = "scenario_id")) %>% 
  mutate(
    algorithm = factor(algorithm, levels = c(order_algorithms))
  ) %>% 
  unnest(correlation_coeffs, keep_empty = TRUE) %>% 
  group_by(id, algorithm, simulation_id) %>% 
  mutate(
    n_sample = unlist(n_sample),
    l = unlist(l),
    n_dim = lengths(variance),
    main_dim = 1:n(),
    n_main_dim = max(main_dim)
  ) %>% 
  ungroup()

# Compute the number of simulations for each scenario
df_scenario_corr %>% 
  distinct(id, simulation_id, n_sample, l, n_dim, n_main_dim) %>% 
  group_by(n_sample, l, n_dim, n_main_dim) %>% 
  summarise(total = n()) %>% 
  View()

df_scenario_corr %>% 
  ggplot(aes(x = l, y = correlation_coeffs, group = algorithm, colour = algorithm)) + 
  geom_point() + 
  geom_jitter()

df_scenario_corr %>% 
  ggplot(aes(x = l, y = -log(1-correlation_coeffs), group = algorithm, colour = algorithm)) + 
  geom_point() + 
  geom_jitter() +
  xlab("\u2113 value") + 
  ylab("Correlation")

df_scenario_corr %>% 
  filter(l > 50) %>% 
  ggplot(aes(x = l, y = -log(1-correlation_coeffs), group = interaction(l, algorithm))) + 
  geom_boxplot(aes(colour = algorithm), width = 25) + 
  stat_summary(
    aes(group = algorithm, colour = algorithm), 
    fun = median,
    geom="line",
    position = position_dodge(width = 25)
  )

fake_plot <- df_scenario_corr %>% 
  filter(!is.na(correlation_coeffs)) %>% 
  group_by(l, algorithm) %>% 
  summarise(
    mean_corr = mean(correlation_coeffs),
    mean_trans_corr = mean(-log(1-correlation_coeffs))
  ) %>% 
  ggplot(aes(x = l, y = mean_corr, colour = algorithm)) +
  geom_line() +
  geom_point() +
  #xlab("\u2113 value") + 
  ylab("Correlation") +
  theme(
    legend.position="bottom",
    legend.title = element_blank(),
    legend.margin = margin(c(0, 0, 0, 0)),
    legend.key.spacing.x = unit(0, "pt")
  ) +
  guides(colour = guide_legend(nrow = 1))

fake_plot

legend <- cowplot::get_legend(fake_plot)

cairo_pdf(
  "/Users/cristianpachon/MEGA/tesis_cristian/mds_for_big_data/to_ADAC/legend_temp.pdf",
  family = "Helvetica",
  width = 8,
  height = 0.3
)
#cairo_pdf("images/legend.pdf", family = "Helvetica", width = 8, height = 8)
grid.newpage()
grid.draw(legend)
dev.off()

cairo_pdf("images/correlation_l_selection.pdf", family = "Helvetica", width = 4, height = 4)
df_scenario_corr %>% 
  filter(!is.na(correlation_coeffs), l<800) %>% 
  group_by(l, algorithm) %>% 
  summarise(
    mean_corr = mean(correlation_coeffs),
    mean_trans_corr = mean(-log(1-correlation_coeffs))
  ) %>% 
  ggplot(aes(x = l, y = mean_corr, colour = algorithm)) +
  geom_line() +
  geom_point() +
  xlab(expression("\u2113 value")) + 
  theme(legend.position = "none") +
  ylab("Correlation") 
dev.off()

#  Analyse time -----
df_scenario_times <- df_scenarios %>% 
  left_join(df_times, by = c("id" = "scenario_id")) %>% 
  mutate(
    l = unlist(l),
    algorithm = factor(algorithm, levels = c(order_algorithms)),
    n_sample = unlist(n_sample),
    n_dim = lengths(variance),
    n_main_dim = map_int(variance, function(x) sum(x != 1))
  )
  
df_scenario_times %>% 
  #filter(algorithm != "reduced_mds") %>% 
  ggplot(aes(x = l, y = elapsed_time, group = algorithm, colour = algorithm)) + 
  geom_line()

df_scenario_times %>% 
  ggplot(aes(x = l, y = log2(elapsed_time), group = algorithm, colour = algorithm)) + 
  geom_point() + 
  geom_jitter()

# Raw data of the next plot
df_scenario_times %>% 
  #filter(l > 50) %>% 
  group_by(l, algorithm) %>% 
  summarise(
    mean_time = mean(elapsed_time)
  ) %>% View()
  
cairo_pdf("images/time_l_selection.pdf", family = "Helvetica", width = 4, height = 4)
df_scenario_times %>% 
  filter(!is.na(elapsed_time), l<800) %>% 
  #filter(l > 50) %>% 
  group_by(l, algorithm) %>% 
  summarise(
    mean_time = mean(elapsed_time)
  ) %>% 
  ggplot(aes(x = l, y = mean_time, colour = algorithm)) +
  geom_line() +
  geom_point() +
  xlab(expression("\u2113 value")) + 
  theme(legend.position = "none") +
  ylab("Elapsed time (sec.)")
dev.off()

View(df_scenario_times)


# Analysis of eigenvalues ----
# According to Wikipedia (and LyD), bias = E(estimator) - parameter
df_scenario_eigen <- df_scenarios %>% 
  left_join(df_eigenvalues, by = c("id" = "scenario_id")) %>% 
  unnest(eigenvalues) %>% 
  group_by(id, algorithm, simulation_id) %>% 
  mutate(
    l = unlist(l),
    algorithm = factor(algorithm, levels = c(order_algorithms)),
    n_sample = unlist(n_sample),
    n_dim = lengths(variance),
    main_dim = 1:n(),
    n_main_dim = map_int(variance, function(x) sum(x != 1))
  ) %>% 
  ungroup()

View(df_scenario_eigen)

# Sanity checks
df_scenario_eigen %>% 
  distinct(id, simulation_id, algorithm, l, n_main_dim, n_dim) %>% 
  group_by(algorithm, l, n_main_dim, n_dim) %>% 
  summarise(
    total = n(),
    min_total = min(total),
    max_total = max(total)
  ) %>% 
  ungroup() %>% 
  distinct(min_total, max_total)

cairo_pdf("images/rmse_l_selection.pdf", family = "Helvetica", width = 4, height = 4)
df_scenario_eigen %>% 
  filter(l > 50, l<800) %>% 
  group_by(l, algorithm) %>% 
  summarise(
    bias = mean(eigenvalues) - 15,
    rmse = sqrt(var(eigenvalues) + bias^2)
  ) %>% 
  ggplot(aes(x = l, y = rmse, colour = algorithm)) +
  geom_line() +
  geom_point() +
  xlab(expression("\u2113 value")) + 
  #xlab("\u2113 value") + 
  theme(legend.position = "none") +
  ylab(TeX("RMSE"))
dev.off()


library(latex2exp)
TeX("RMSE $\\sigma$")
