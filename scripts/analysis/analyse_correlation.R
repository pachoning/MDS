library(tidyverse)
library(ggplot2)

# Load data ----
data_path = file.path(getwd(), 'data')
load(file.path(data_path, "df_scenarios_full.RData"))
load(file.path(data_path, "df_correlation_full.RData"))

# Manipulate data ----
scenario_identifier = c("sample_size", "n_cols", "n_main_dimensions", "sd_main")

# Avoid using scenarions which sample size is 10^6
scenarios_with_main_dimesions = df_scenarios_full %>% filter(n_main_dimensions > 0, !is.na(processed_at))

# Join scenarios and correlation
df_join_scenarios_correlation = scenarios_with_main_dimesions %>% 
  select_at(c("id", scenario_identifier))   %>% 
  left_join(
    df_correlation_full %>% select(-n_main_dimensions),
    by = c("id" = "scenario_id")
  ) %>% unnest(correlation_vector) %>% 
  rename(correlation = correlation_vector) %>% 
  group_by(id, method_name, num_sim) %>% 
  mutate(
    dim = paste0("dim_", sprintf("%02d", 1:n())),
    log_1_correlation = log(1-correlation)
  )

# Analyse data  ----
# Plots
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
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          legend.position="none")
  
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

# Plot of the correlation of the first dimension
df_join_scenarios_correlation %>% 
  group_by(method_name, sample_size, dim) %>%
  summarise(mean_correlation = mean(correlation)) %>%
  mutate(sample_size = as.factor(sample_size)) %>% 
  ggplot(aes(x = sample_size, y = mean_correlation, color = method_name, group = method_name)) +
  geom_point() +
  facet_wrap( ~ dim, ncol = 2) +
  theme(panel.spacing.y=unit(0.5, "lines")) +
  ggsave("/Users/cristianpachongarcia/Documents/phd/papers/mds_for_big_data/images/correlation_all.png", 
         dpi=300, dev='png', height=8, width=10, units="in")

df_join_scenarios_correlation %>% 
  filter(id %in% scenarios_with_main_dimesions$id, sample_size<10^5) %>% 
  group_by(method_name, sample_size, dim) %>%
  summarise(mean_correlation = mean(correlation)) %>%
  mutate(sample_size = as.factor(sample_size)) %>% 
  ggplot(aes(x = sample_size, y = mean_correlation, color = method_name, group = method_name)) +
  geom_point() +
  facet_wrap( ~ dim, ncol = 2) +
  theme(panel.spacing.y=unit(0.5, "lines")) +
  ggsave("/Users/cristianpachongarcia/Documents/phd/papers/mds_for_big_data/images/correlation_small.png", 
         dpi=300, dev='png', height=8, width=6.5, units="in")

df_join_scenarios_correlation %>% 
  filter(id %in% scenarios_with_main_dimesions$id, sample_size>=10^5) %>% 
  group_by(method_name, sample_size, dim) %>%
  summarise(mean_correlation = mean(correlation)) %>%
  mutate(sample_size = as.factor(sample_size)) %>% 
  ggplot(aes(x = sample_size, y = mean_correlation, color = method_name, group = method_name)) +
  geom_point() +
  facet_wrap( ~ dim, ncol = 2) +
  theme(panel.spacing.y=unit(0.5, "lines")) +
  ggsave("/Users/cristianpachongarcia/Documents/phd/papers/mds_for_big_data/images/correlation_big.png", 
         dpi=300, dev='png', height=8, width=6.5, units="in")


# Particular cases
df_join_scenarios_correlation %>% 
  filter(sample_size == 10^6, n_main_dimensions == 10, n_cols == 100) %>% 
  ggplot(aes(x = method_name, y = correlation, fill = method_name)) + 
  geom_boxplot() + 
  facet_wrap(. ~ dim, ncol = 5) +
  theme(
    panel.spacing.y=unit(0.5, "lines"), 
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  ) +
  ggsave("/Users/cristianpachongarcia/Documents/phd/papers/mds_for_big_data/images/correlation_1000000_100_10.png", 
         dpi=300, dev='png', height=8, width=10, units="in")
