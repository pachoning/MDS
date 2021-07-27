library(bigmds)
library(tidyverse)
library(RColorBrewer)

# Load data ----
load("data/EMNIST/divide.RData")
load("data/EMNIST/fast.RData")
load("data/EMNIST/gower.RData")
load("data/EMNIST/all_data.RData")

# Variables ----
target_filter <- c("0", "1", "a", "b", "A", "B")
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

# Analysis ----
df_aux <- df_mds_fast
df_aux$type_data <- type_data

df_aux %>% 
  mutate(category = ifelse(
    target %in% as.character(0:9),
    "Number",
    ifelse(
      target %in% letters,
      "Small_letter",
      ifelse(
        target %in% toupper(letters),
        "Upper_letter",
        "UNK"
      )
    )
  )) %>% 
  group_by(category, type_data) %>% 
  summarise(n())

df_aux %>% 
  group_by(type_data) %>% 
  summarise(n())

# Correlation and mean of MDS ----
df_corr <- df_mds_fast
cor(df_corr[, c(1,2)])
apply(df_corr[, c(1,2)], MARGIN = 2,mean)


# Sampling ----
idxs <- 1:nrow(df_mds_divide)
target <- df_mds_divide$target
idxs_filtered <- idxs[target %in% target_filter]
idx_sample <- sample(idxs_filtered, size = sample_size)
idx_sample <- sort(idx_sample)

df_mds_divide_sample <- df_mds_divide[idx_sample, ]
df_mds_fast_sample <- df_mds_fast[idx_sample, ]
df_mds_gower_sample <- df_mds_gower[idx_sample, ]


# Some validations ----
# Check the correlation and the mean
cor(df_mds_divide[, c(1,2), drop = FALSE])
apply(df_mds_divide[, c(1,2), drop = FALSE], MARGIN = 2, mean)

cor(df_mds_fast[, c(1,2), drop = FALSE])
apply(df_mds_fast[, c(1,2), drop = FALSE], MARGIN = 2, mean)

cor(df_mds_gower[, c(1,2), drop = FALSE])
apply(df_mds_gower[, c(1,2), drop = FALSE], MARGIN = 2, mean)

# Plots ----
color_fig <- c(brewer.pal(5, "Set1"), "black")

# Embeddings
df_mds_divide_sample %>% 
  ggplot(aes(x = x, y = y, color = target, tit)) + 
  geom_point() + 
  scale_color_manual(values =color_fig) +
  ggsave(file.path(getwd(), "images", "emnist_divide.png"))

df_mds_gower_sample %>% 
  ggplot(aes(x = x, y = -y, color = target, tit)) + 
  geom_point() +
  scale_color_manual(values = color_fig) +
  ggsave(file.path(getwd(), "images", "emnist_gower.png"))

df_mds_fast_sample %>% 
  ggplot(aes(x = -x, y = y, color = target, tit)) + 
  geom_point() +
  scale_color_manual(values = color_fig) +
  ggsave(file.path(getwd(), "images", "emnist_fast.png"))

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
