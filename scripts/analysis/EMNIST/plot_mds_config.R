library(bigmds)
library(tidyverse)

# Load data ----
load("data/EMNIST/divide.RData")
load("data/EMNIST/fast.RData")
load("data/EMNIST/gower.RData")

# Variables ----
target_filter <- as.character(0:4)
sample_size <- NULL
unique_categories <- df_mds_divide %>% count(target)

images_per_category <- df_mds_divide %>% 
  group_by(target) %>% 
  summarise(total = n())
  
if (is.null(target_filter)) {
  target_filter <- unique(df_mds_divide$target)
}

if (is.null(sample_size)) {
  sample_size <- images_per_category %>%
    filter(target %in% target_filter) %>%
    summarise(total = sum(total)) %>% 
    pull()
}

# Sampling ----
idxs <- 1:nrow(df_mds_divide)
target <- df_mds_divide$target
idxs_filtered <- idxs[target %in% target_filter]
idx_sample <- sample(idxs_filtered, size = sample_size)
idx_sample <- sort(idx_sample)

df_mds_divide_sample <- df_mds_divide[idx_sample, ]
df_mds_fast_sample <- df_mds_fast[idx_sample, ]
df_mds_gower_sample <- df_mds_gower[idx_sample, ]


# Plots ----
# Embeddings
df_mds_divide_sample %>% 
  ggplot(aes(x = x, y = y, color = target, tit)) + 
  geom_point() +
  ggtitle("Divide")

df_mds_fast_sample %>% 
  ggplot(aes(x = x, y = y, color = target, tit)) + 
  geom_point() +
  ggtitle("Fast")

df_mds_gower_sample %>% 
  ggplot(aes(x = x, y = y, color = target, tit)) + 
  geom_point() +
  ggtitle("Gower")

# Correlations
cor(df_mds_divide_sample$x, df_mds_fast_sample$x)
cor(df_mds_divide_sample$x, df_mds_gower_sample$x)
cor(df_mds_fast_sample$x, df_mds_gower_sample$x)

cor(df_mds_divide_sample$y, df_mds_fast_sample$y)
cor(df_mds_divide_sample$y, df_mds_gower_sample$y)
cor(df_mds_fast_sample$y, df_mds_gower_sample$y)

# Time ----
elapsed_time_divide
elapsed_time_gower
elapsed_time_fast
