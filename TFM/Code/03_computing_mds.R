################################################################################
# The goal of this script is to do MDS using the data splitted
################################################################################


# 01 Load libraries ----
library(gower)
library(cluster)
library(tidyverse)

# 02 Initial parameters ----
input_directory = "data/Bike-Sharing-Dataset/"
input_file = "df_split.RData"
output_directory = ""



# 03 Load the data ----
load(paste0(input_directory, input_file))
str(df_split)

# 03 Compute Gower distance ----
# We don't consider the following variables neither to compute the distances
# nor to fo the MDS:
#   instant
#   dteday

df_cols_selection = df_split %>% 
  select(
    - instant,
    - dteday
  )


gower_distance = df_cols_selection %>% 
  daisy(
    metric = "gower"
  ) %>% 
  as.matrix()



# 04 MDS ----
# k is the number of dim
fit <- cmdscale(
  gower_distance, 
  eig=TRUE, 
  k=2
) 
fit


# 05 Plot data ----
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(
  x, 
  y, 
  xlab="Coordinate 1", 
  ylab="Coordinate 2", 
  main="Metric MDS",	
  type="n"
)
text(x, y, labels = row.names(df_total), cex=.7)
