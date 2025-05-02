#library(bigmds)
library(tidyverse)
library(RColorBrewer)

# Load data ----
load("data/EMNIST/interpolation.RData")
load("data/EMNIST/divide.RData")
load("data/EMNIST/fast.RData")
load("data/EMNIST/landmark.RData")
load("data/EMNIST/pivot.RData")
load("data/EMNIST/reduced.RData")
load("data/EMNIST/all_data.RData")

# Variables ----
target_filter <- c("0", "1", "S", "r")
#target_filter <- c("0", "1", "a", "b", "A", "B")
sample_size <- 2000
unique_categories <- df_mds_divide %>% count(target)

images_per_category <- df_mds_divide %>% 
  group_by(target) %>% 
  summarise(
    total = n()
  )

images_per_category$rate <- images_per_category$total/sum(images_per_category$total)
 
images_per_category %>%
  filter(target %in% target_filter) %>% 
  View()

if (is.null(target_filter)) {
  target_filter <- unique(df_mds_divide$target)
}

if (is.null(sample_size)) {
  sample_size <- images_per_category %>%
    filter(target %in% target_filter) %>%
    summarise(total = sum(total)) %>% 
    pull()
}

if (FALSE){
  
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
}

# Sampling ----
idxs <- 1:nrow(df_mds_divide)
target <- df_mds_divide$target
idxs_filtered <- idxs[target %in% target_filter]
idx_sample <- sample(idxs_filtered, size = sample_size)
idx_sample <- sort(idx_sample)

df_mds_interpolation_sample <- df_mds_interpolation[idx_sample, ]
df_mds_divide_sample <- df_mds_divide[idx_sample, ]
df_mds_fast_sample <- df_mds_fast[idx_sample, ]
df_mds_landmark_sample <- df_mds_landmark[idx_sample, ]
df_mds_pivot_sample <- df_mds_pivot[idx_sample, ]
df_mds_reduced_sample <- df_mds_reduced[idx_sample, ]


if (FALSE) {
  
  # Some validations ----
  # Check the correlation and the mean
  sum(apply(all_data_pixels, 2, var))
  
  cor(df_mds_interpolation[, c(1,2), drop = FALSE])
  cov(df_mds_interpolation[, c(1,2), drop = FALSE])
  cov(df_mds_interpolation[1:150, c(1,2), drop = FALSE])
  sum(diag(cov(df_mds_interpolation[, c(1,2)])))
  apply(df_mds_interpolation[, c(1,2), drop = FALSE], MARGIN = 2, mean)
  
  cor(df_mds_divide[, c(1,2), drop = FALSE])
  apply(df_mds_divide[, c(1,2), drop = FALSE], MARGIN = 2, mean)
  
  cor(df_mds_fast[, c(1,2), drop = FALSE])
  apply(df_mds_fast[, c(1,2), drop = FALSE], MARGIN = 2, mean)
  
  cor(df_mds_landmark[, c(1,2), drop = FALSE])
  cov(df_mds_landmark[, c(1,2), drop = FALSE])
  apply(df_mds_landmark[, c(1,2), drop = FALSE], MARGIN = 2, mean)
  
  cor(df_mds_pivot[, c(1,2), drop = FALSE])
  apply(df_mds_pivot[, c(1,2), drop = FALSE], MARGIN = 2, mean)
  
  cor(df_mds_reduced[, c(1,2), drop = FALSE])
  apply(df_mds_reduced[, c(1,2), drop = FALSE], MARGIN = 2, mean)
  

  # Embeddings
  df_mds_interpolation_sample %>% 
    ggplot(aes(x = x, y = -y, color = target, tit)) + 
    geom_point() +
    scale_color_manual(values = color_fig) +
    ggsave(file.path(getwd(), "images", "emnist_gower_new.png"))
  
  df_mds_divide_sample %>% 
    ggplot(aes(x = x, y = y, color = target, tit)) + 
    geom_point() + 
    scale_color_manual(values =color_fig) +
    ggsave(file.path(getwd(), "images", "emnist_divide_new.png"))
  
  df_mds_fast_sample %>% 
    ggplot(aes(x = -x, y = y, color = target, tit)) + 
    geom_point() +
    scale_color_manual(values = color_fig) +
    ggsave(file.path(getwd(), "images", "emnist_fast_new.png"))
  
  df_mds_landmark_sample %>% 
    ggplot(aes(x = -x, y = y, color = target, tit)) + 
    geom_point() +
    scale_color_manual(values = color_fig) +
    ggsave(file.path(getwd(), "images", "emnist_landmark_new.png"))
  
  df_mds_pivot_sample %>% 
    ggplot(aes(x = -x, y = y, color = target, tit)) + 
    geom_point() +
    scale_color_manual(values = color_fig) +
    ggsave(file.path(getwd(), "images", "emnist_pivot_new.png"))
  
  df_mds_reduced_sample %>% 
    ggplot(aes(x = -x, y = y, color = target, tit)) + 
    geom_point() +
    scale_color_manual(values = color_fig) +
    ggsave(file.path(getwd(), "images", "emnist_reduced_new.png"))
  
  # Correlations
  cor(df_mds_divide_sample$x, df_mds_fast_sample$x)
  cor(df_mds_divide_sample$x, df_mds_interpolation_sample$x)
  cor(df_mds_fast_sample$x, df_mds_interpolation_sample$x)
  
  cor(df_mds_divide_sample$y, df_mds_fast_sample$y)
  cor(df_mds_divide_sample$y, df_mds_interpolation_sample$y)
  cor(df_mds_fast_sample$y, df_mds_interpolation_sample$y)
}
# Plots ----
color_fig <- c(brewer.pal(5, "Set1"), "black")

# Embedding all
df_mds_interpolation_sample$algorithm <- "interpolation_mds"
df_mds_divide_sample$algorithm <- "divide_conquer_mds"
df_mds_fast_sample$algorithm <- "fast_mds"
df_mds_landmark_sample$algorithm <- "landmark_mds"
df_mds_pivot_sample$algorithm <- "pivot_mds"
df_mds_reduced_sample$algorithm <- "reduced_mds"

df_mds_landmark_sample$y <- -df_mds_landmark_sample$y
#df_mds_interpolation_sample$x <- -df_mds_interpolation_sample$x
#df_mds_reduced_sample$y <- -df_mds_reduced_sample$y
df_mds_pivot_sample$y <- -df_mds_pivot_sample$y
#df_mds_divide_sample$x <- -df_mds_divide_sample$x
df_mds_fast_sample$x <- -df_mds_fast_sample$x
df_mds_fast_sample$y <- -df_mds_fast_sample$y



df_embedding <- rbind(
  df_mds_interpolation_sample,
  df_mds_divide_sample,
  df_mds_fast_sample,
  df_mds_landmark_sample,
  df_mds_pivot_sample,
  df_mds_reduced_sample
)

order_algorithms <- c(
  "landmark_mds",
  "interpolation_mds",
  "reduced_mds",
  "pivot_mds",
  "divide_conquer_mds",
  "fast_mds"
)

df_embedding$algorithm <- factor(df_embedding$algorithm, levels = c(order_algorithms))

pdf('images/emnist_all.pdf', width = 8, height = 4)
#png('images/emnist_all.png')
df_embedding %>% 
  ggplot(aes(x = x, y = y, color = target, tit)) + 
  geom_point(size = .8) +
  facet_wrap(vars(algorithm)) +
  scale_color_manual(
    values = c("black", "blue", "red", "green3")
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.background = element_rect(fill='transparent', color=NA)
  )
dev.off()


# Time ----
elapsed_time_landmark
elapsed_time_interpolation
elapsed_time_reduced
elapsed_time_pivot
elapsed_time_divide
elapsed_time_fast

eigen_landmark
eigen_interpolation
eigen_reduced
eigen_pivot
eigen_divide
eigen_fast
