library(tidyverse)
library(ggplot2)
library(ggridges)

# Load data  ----
data_path = file.path(getwd(), 'data')
load(file.path(data_path, "df_scenarios_full.RData"))
load(file.path(data_path, "df_time_full.RData"))

# Manipulate data ----
scenario_identifier = c("sample_size", "n_cols", "n_main_dimensions", "var_main")

# Avoid using scenarions which sample size is 10^6
df_scenarios_full_filtered = df_scenarios_full %>% filter(!is.na(processed_at))

# Join scenarios and time
df_join_scenarios_time = df_scenarios_full_filtered %>% 
  select_at(c("id", scenario_identifier)) %>% 
  left_join(
    df_time_full,
    by = c("id" = "scenario_id")
  ) %>% 
  mutate(n_main_dimensions=as.factor(n_main_dimensions),
         n_cols = as.factor(n_cols),
         log_elapsed_time = log(elapsed_time))

# Analyse data  ----
# Anova
df_anova = df_join_scenarios_time %>% mutate(sample_size = as.factor(sample_size))
linear_model = lm(elapsed_time ~ n_cols + sample_size + n_main_dimensions + method_name, data = df_anova)
anova(linear_model)
summary(linear_model)

# Print a density plot and save it for each scenario
total_scenario = nrow(df_scenarios_full_filtered)
for(i_scenario in 1:total_scenario){
  current_scenario = df_scenarios_full_filtered[i_scenario, ]
  results = df_join_scenarios_time[df_join_scenarios_time$id == current_scenario$id, ]
  
  plot_title = paste0("sample_size: ", current_scenario$sample_size, "; n_cols: ", current_scenario$n_cols, 
                      "; n_main_dim: ", current_scenario$n_main_dimensions, "; var_main: ", current_scenario$var_main)
  plot_name = paste0("sample_size:", current_scenario$sample_size, "__n_cols:", current_scenario$n_cols, 
                     "; __n_main_dim:", current_scenario$n_main_dimensions, 
                     "; __var_main:",current_scenario$var_main, ".png")
  
  p = ggplot(results, aes(elapsed_time, group = method_name, color = method_name)) +
    geom_density() + 
    ggtitle(plot_title) +
    theme(plot.title = element_text(hjust = 0.5))
  
  png(file.path("plots", "time", plot_name))
  print(p)
  dev.off()
}

# Mean/Median/... for each scenario
main_statistics = df_join_scenarios_time %>% 
  group_by_at(c(scenario_identifier, "method_name")) %>% 
  summarise(
    mean_time = mean(elapsed_time)
  ) %>% 
  pivot_wider(names_from = method_name, values_from = mean_time)

methods_names = colnames(main_statistics)[-4:-1]
main_statistics$best_method = methods_names[apply(main_statistics[,methods_names],1,which.min)]
main_statistics %>% View

# Pots by sample size
df_summary_sample_size = df_join_scenarios_time %>%
  group_by(sample_size, method_name) %>%
  summarise(
    mean_elapsed_time = mean(elapsed_time), 
    log_mean_elapsed_time = log(mean_elapsed_time),
    mean_log_elapsed_time = mean(log_elapsed_time)
  )

df_summary_sample_size %>% 
  mutate(sample_size = as.factor(sample_size)) %>% 
  ggplot(aes(x = sample_size, y = mean_log_elapsed_time, group = method_name, color = method_name)) +
  geom_point(size = 2) +
  geom_line() +
  theme(panel.spacing.y=unit(0.5, "lines")) +
  ggsave("~/Documents/phd/papers/mds_for_big_data/images/mean_log_time_all.png", 
         dpi=300, dev='png', height=8, width=10, units="in")


df_summary_sample_size %>% 
  filter(sample_size < 10^5)%>% 
  ggplot(aes(x = sample_size, y = mean_elapsed_time, group = method_name, color = method_name)) +
  geom_point(size = 2) +
  geom_line(alpha=0.5) +
  theme(
    panel.spacing.y=unit(0.5, "lines"), 
    axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  ) +
  ggsave("~/Documents/phd/papers/mds_for_big_data/images/time_small.png", 
         dpi=300, dev='png', height=8, width=6.5, units="in")


df_summary_sample_size %>% 
  filter(sample_size >= 10^5)%>% 
  ggplot(aes(x = sample_size, y = mean_elapsed_time, group = method_name, color = method_name)) +
  geom_point(size = 2) +
  geom_line(alpha=0.5) +
  theme(
    panel.spacing.y=unit(0.5, "lines"), 
    axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  ) +
  ggsave("~/Documents/phd/papers/mds_for_big_data/images/time_big.png", 
         dpi=300, dev='png', height=8, width=6.5, units="in")


# Time for a particular case
df_join_scenarios_time %>% 
  filter(sample_size == 10^6, n_main_dimensions == 10, n_cols == 100) %>% 
  group_by(method_name) %>% 
  summarise(p1 = quantile(elapsed_time, probs = 0.05/2),
            p2 = quantile(elapsed_time, probs = 1-0.05/2))

