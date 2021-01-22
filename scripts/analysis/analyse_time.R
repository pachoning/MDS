library(tidyverse)
library(ggplot2)

# Load data  ----
data_path = file.path(getwd(), 'data')
load(file.path(data_path, "df_scenarios_full.RData"))
load(file.path(data_path, "df_time_full.RData"))

# Manipulate data ----
scenario_identifier = c("sample_size", "n_cols", "n_main_dimensions", "sd_main")

# Avoid using scenarions which sample size is 10^6
df_scenarios_full_filtered = df_scenarios_full %>% 
  filter(sample_size != 10^6)

# Join scenarios and time
df_join_scenarios_time = df_scenarios_full_filtered %>% 
  select_at(c("id", scenario_identifier)) %>% 
  left_join(
    df_time_full,
    by = c("id" = "scenario_id")
  ) %>% 
  mutate(n_main_dimensions=as.factor(n_main_dimensions), 
         sample_size=as.factor(sample_size),
         n_cols = as.factor(n_cols),
         log_elapsed_time = log(elapsed_time))

# Analyse data  ----
# Anova
linear_model = lm(log_elapsed_time ~ n_cols + sample_size + n_main_dimensions + method_name, data = df_join_scenarios_time)
anova(linear_model)
summary(linear_model)

# Print a density plot and save it for each scenario
total_scenario = nrow(df_scenarios_full_filtered)
for(i_scenario in 1:total_scenario){
  current_scenario = df_scenarios_full_filtered[i_scenario, ]
  results = df_join_scenarios_time[df_join_scenarios_time$id == current_scenario$id, ]
  
  plot_title = paste0("sample_size: ", current_scenario$sample_size, "; n_cols: ", current_scenario$n_cols, 
                      "; n_main_dim: ", current_scenario$n_main_dimensions, "; sd_main: ", current_scenario$sd_main)
  plot_name = paste0("sample_size:", current_scenario$sample_size, "__n_cols:", current_scenario$n_cols, 
                     "; __n_main_dim:", current_scenario$n_main_dimensions, 
                     "; __sd_main:",current_scenario$sd_main, ".png")
  
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
main_statistics %>% View
