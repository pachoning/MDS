################################################################################
# The goal of this script is to split the data/create a column that indicates
# the set an observation belogs. The output is the same data frame with an 
# additional column indicating the group
################################################################################

#01 Load libraries and other sources ----
library(lubridate)
library(tidyverse)
library(gtools)

source("Tools/R/split_sample.R")

# 02 Initial parameters ----
input_directory = file.path(
  getwd(),
  "Data",
  "Bike-Sharing-Dataset"
)
input_file = "df_input.RData"


output_directory = file.path(
  getwd(),
  "Data",
  "Bike-Sharing-Dataset"
)
  


# 03 Load data ----
load(file.path(input_directory, input_file))
str(df_input)

# 04 Split randomnly----
df_split = split_sample(
  x = df_input,
  is_random_method = TRUE,
  number_groups = 1000,
  variable = NULL,
  return_data_frame = TRUE
)

df_split_random = df_split


# 05 Split as time series ----
df_split = split_sample(
  x = df_input,
  is_random_method = FALSE,
  number_groups = 5,
  variable = "dteday",
  return_data_frame = TRUE
)

df_split %>% 
  group_by(
    group_member
  ) %>% 
  summarise(
    total_obs = n(),
    min_value = min(dteday),
    max_value = max(dteday)
  )


# 06 Split using season variable ----

df_split = split_sample(
  x = df_input,
  is_random_method = FALSE,
  number_groups = NULL,
  variable = "season",
  return_data_frame = TRUE
)
df_split %>% 
  group_by(
    group_member
  ) %>% 
  summarise(
    total_obs = n(),
    min_value = min(season),
    max_value = max(season)
  )


# 07 Split using a continuous variable ----
df_split = split_sample(
  x = df_input,
  is_random_method = FALSE,
  number_groups = 5,
  variable = "windspeed",
  return_data_frame = TRUE
)
df_split %>% 
  group_by(
    group_member
  ) %>% 
  summarise(
    total_obs = n(),
    min_value = min(windspeed),
    max_value = max(windspeed)
  )

# 08 save ----
df_split = df_split_random
save(
  df_split,
  file = file.path(output_directory, "df_split.RData")
)

