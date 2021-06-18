library(tidyverse)
load("data/cmdscale/df_cost_time.RData")

df_cost_time %>%
  ggplot(aes(x = l, y = time, color = event)) +
  geom_line()

df_cost_time %>% filter(l == 2500)
