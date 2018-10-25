# 01 Load libraries and other sources ----
library(lubridate)
library(tidyverse)
library(gtools)

# 02 variables ----
path_input_folder = file.path(
  getwd(), 
  "Data",
  "Bike-Sharing-DataSet"
)

input_file_name = "df_input.RData"


path_input_file = file.path(
  path_input_folder,
  input_file_name
  
)

load(path_input_file)

df_input  


# 03 Some plots ----
# Period of time recorded
range(df_input$dteday)
range(df_input$instant)
range(df_input$casual)
range(df_input$registered)
unique(df_input$hr)


# users
ggplot() + 
  geom_density(aes(x = casual), colour="red", data = df_input) + 
  geom_density(aes(x = registered), colour = "blue", data = df_input)


# Time series record per hour
df_input %>% 
  mutate(
    total_users = casual + registered
  ) %>% 
  ggplot(
    aes(
      x = hr,
      y = total_users
    )
  ) +
  geom_boxplot() 
 


# weathersit
df_input %>% 
  mutate(
    total_users = casual + registered
  ) %>% 
  ggplot(
    aes(
      x = weathersit,
      y = total_users
    )
  ) +
  geom_boxplot()

