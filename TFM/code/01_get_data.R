################################################################################
# The goal of this script is to get the input data and generate a RData as 
# output with the correct format for each column
################################################################################

# 01 Load libraries and other sources ----
library(lubridate)
library(tidyverse)
library(gtools)


# 02 Parameters/Variables ----
input_directory = "Data/Bike-Sharing-Dataset/"
input_file = "hour.csv"
output_directory = "data/Bike-Sharing-Dataset/"

# 03 input data ----
df = read_csv(paste0(input_directory, input_file))



# 04 Manipulate data ----
df_input = df %>%
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


# 05 Save ----
save(
  df_input,
  file = paste0(output_directory, "df_input.RData")
)
