df_join = df_time %>% 
  left_join(
    df_scenarios %>% select(id, n_main_dimensions),
    by = c("scenario_id" = "id")
  )

df_join %>% 
  mutate(n_main_dimensions=as.factor(n_main_dimensions)) %>% 
  group_by(n_main_dimensions, method_name) %>% 
  summarise(mean_time = mean(elapsed_time))%>% 
  ggplot(aes(x=n_main_dimensions, y=mean_time, color=method_name, group=method_name)) +
  geom_line()


df_scenarios
df_time %>% View
