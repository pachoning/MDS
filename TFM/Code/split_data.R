#01 Load libraries and other sources ----
library(lubridate)
library(tidyverse)
library(gtools)


source("Tools/R/split_sample.R")


#02 input data ----
path_file = "Data/Bike-Sharing-Dataset/hour.csv"
df = read_csv(path_file)


#03 Manipulate data ----
df = df %>%
  mutate(
    dteday = lubridate::date(dteday),
    season = as.ordered(season),
    yr = as.ordered(yr),
    mnth = as.ordered(mnth),
    hr = as.ordered(hr),
    holiday = as.factor(holiday),
    weekday = as.ordered(weekday),
    workingday = as.factor(workingday),
    weathersit = as.ordered(weathersit)
  )

#04 Split ----

# Randomly
df_split = split_sample(
  x = df,
  is_random_method = TRUE,
  number_groups = 5,
  variable = NULL,
  return_data_frame = TRUE
)



# Time series
df_split = split_sample(
  x = df,
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


# Season variable

df_split = split_sample(
  x = df,
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


# Continuous variable
df_split = split_sample(
  x = df,
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
