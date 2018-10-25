################################################################################
# The goal of this script is to do MDS using the data splitted
################################################################################

# 01 Load libraries ----
library(tidyverse)
source("tools/divide_conquer_mds.R")
library(ggplot2)
library(ggrepel)


# 02 Testing the algorithm on Australian cities ----
# Traditional MDS
url <- "http://rosetta.reltech.org/TC/v15/Mapping/data/dist-Aus.csv"
dist.au <- read.csv(url)
dist.au

row.names(dist.au) <- dist.au[, 1]
dist.au <- dist.au[, -1]
dist.au

fit <- cmdscale(dist.au, eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]

plot(
  x, 
  y, 
  main = "using traditional mds",
  xlab = "",
  ylab = "",
  pch = 9, 
  xlim = range(x) + c(0, 600) 
)


city.names <- c(
  "Adelaide", "Alice Springs", "Brisbane", "Darwin", "Hobart", 
  "Melbourne", "Perth", "Sydney"
)

text(x, y, pos = 4, labels = city.names)

# Divide and conquer algorithm
fit2 <- cmdscale(dist.au, eig = TRUE, k = 4)
dist.au
x_corrdinates <- fit2$points

# This should be similar to dist.au
daisy(
  x_corrdinates,
  metric = "euclidean"
)

dist.au

groups = c(1,1,1,1,2,2,2,2)
number_coordinates = 2

mds_algorithm = divide_conquer_mds(
  x = x_corrdinates,
  groups = groups,
  number_coordinates = 2,
  metric = "euclidean"
)

plot(
  mds_algorithm$mds[,1], 
  mds_algorithm$mds[,2], 
  pch = 19, 
  xlim = range(x) + c(0, 600),
  main = "using divide-conquer-mds",
  xlab = "",
  ylab = ""
)

city.names <- c(
  "Adelaide", "Alice Springs", "Brisbane", "Darwin", "Hobart", 
  "Melbourne", "Perth", "Sydney"
)

text(
  mds_algorithm$mds[,1], 
  mds_algorithm$mds[,2], 
  pos = 4, 
  labels = city.names
)

mds_algorithm$error



# 03 Testing the algorithm on Eurpean cities ----
# Traditional MDS
mds_europe <- cmdscale(eurodist)
# plot(mds_europe, type = 'n')
# text(mds_europe[, 1], mds_europe[, 2], labels(eurodist))

ggplot(
  as.data.frame(mds_europe), 
  aes(
    V1, 
    -V2, 
    label = rownames(mds_europe)
  )
) +
  geom_point() + 
  geom_text_repel() + 
  xlab('x') + 
  ylab('y') +
  labs(
    title = 'Traditional MDS',
    x = 'x',
    y = 'y'
  )

# Divide and conquer algorithm
fit2 <- mds_europe

x_corrdinates <- fit2
dim(x_corrdinates)
# This should be similar to dist.au
daisy(
  x_corrdinates,
  metric = "euclidean"
)
eurodist

groups = c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2)
number_coordinates = 2

mds_algorithm = divide_conquer_mds(
  x = x_corrdinates,
  groups = groups,
  number_coordinates = 2,
  metric = "euclidean"
)

mds_algorithm$mds
mds_algorithm$error

ggplot(
  as.data.frame(mds_algorithm$mds), 
  aes(
    V1, 
    -V2, 
    label = rownames(mds_europe)
  )
) +
  geom_point() + 
  geom_text_repel() + 
  xlab('x') + 
  ylab('y') +
  labs(
    title = 'Divide-conquer-MDS with 2 groups',
    x = 'x',
    y = 'y'
  )

# 04 Applying to a big data set ----
input_directory = file.path(getwd(),"Data", "Bike-Sharing-Dataset")
input_file = "df_split.RData"


load(file.path(input_directory, input_file))
df_split

x_selected = df_split %>% 
  select(
    -instant,
    -dteday,
    -group_member
  )

x_groups = df_split$group_member

mds_algorithm = divide_conquer_mds(
  x = x_selected,
  groups = x_groups,
  number_coordinates = 2,
  metric = "gower"
)

dim(mds_algorithm$mds)
mds_algorithm$error

