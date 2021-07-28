library(tidyverse)

# Load data ----
load("data/EMNIST/df_emnist_gof.RData")
load("data/experiments_processed/df_conversion.RData")

# Manipulate data ----
df_emnist_gof[56, "GOF"] <- 0.73
df_emnist_gof[68, "GOF"] <- 0.78
df_emnist_gof[74, "GOF"] <- 0.79
df_emnist_gof[80,  "GOF"] <- 0.81

df_join <- df_emnist_gof %>% 
  left_join(df_conversion)
levels(df_join$Algorithm) <- c("D&C", "Interp", "Fast")

# Plots ----
df_join %>% 
  ggplot(aes(x = k, y = 100*GOF, color = Algorithm, group = Algorithm)) +
  geom_point(size = 2) +
  geom_line() +
  geom_hline(yintercept = 80, linetype = "dashed") +
  theme(panel.spacing.y=unit(0.5, "lines"), legend.position="bottom") +
  xlab("Dominant dimension") + 
  ylab("GOF") +
  scale_color_manual(values = c("#0000FF", "#FF0000", "#00AF91")) +
  ggsave(file.path(getwd(), "images", "gof_emnist.png"), dev = 'png', width = 22, height = 18, units = "cm")
