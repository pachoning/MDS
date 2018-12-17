################################################################################
# The goal of this script is to do MDS using the data splitted
################################################################################

# 01 Load libraries ----
source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")
set.seed(12345)

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



rotation_matrix = matrix(
  c(cos(180), sin(180), -sin(180), cos(180)), 
  byrow = TRUE,
  nrow = 2
)
  
mds_europe_rotated =   mds_europe %*% rotation_matrix
ggplot(
  as.data.frame(mds_europe_rotated), 
  aes(
    V1, 
    -V2, 
    label = rownames(mds_europe_rotated)
  )
) +
  geom_point() + 
  geom_text_repel() + 
  xlab('x') + 
  ylab('y') +
  labs(
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

# 05 Applying to random data ----
n_obs = 3*10^3
x = data.frame(
  x1 = rnorm(n_obs),
  x2 = rnorm(n_obs),
  x3 = rnorm(n_obs),
  x4 = rnorm(n_obs),
  x5 = rnorm(n_obs),
  x6 = rnorm(n_obs)
)

groups = sample(x = 5, size = nrow(x), replace = TRUE)

mds_algorithm = divide_conquer_mds(
  x = x,
  groups = groups,
  number_coordinates = 2,
  metric = "euclidean"
)


distance_x = daisy(
  x = x,
  metric = "euclidean"
)


mds_x = cmdscale(
  d = distance_x, 
  eig = TRUE, 
  k = 2
)

mds_x = mds_x$points
head(mds_algorithm$mds, 10)
head(mds_x, 10)


rotate_matrix_divide = pracma::procrustes(
  mds_x,
  mds_algorithm$mds
)

mds_rotated = rotate_matrix_divide$P

head(mds_rotated)
head(mds_x)

df_divide_conquer = as.data.frame(mds_rotated)

df_mds_x = as.data.frame(mds_x)
df_divide_conquer$source = "divide"
df_mds_x$source = "classical"

df_divide_conquer$obs = row.names(df_divide_conquer)
df_mds_x$obs = row.names(df_mds_x)

all_results = rbind(
  df_divide_conquer[1:10, ],
  df_mds_x[1:10,]
)

all_results %>% 
  ggplot(aes(x = V1, y = V2, group = source, color = source)) +
  geom_point() +
  geom_text(aes(label=obs),hjust=0, vjust=0)


# 06 Applying to a big data set ----
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

