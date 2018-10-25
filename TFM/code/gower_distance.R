# 00 Packages ----
# Gower package
# https://cran.r-project.org/web/packages/gower/vignettes/intro.html
if(require(gower) == FALSE){
  install.packages("gower")
  library(gower)
} 

# Cluster package
if(require(cluster) == FALSE){
  install.packages("cluster")
  library(cluster)
}


# 01 Getting data ----
dat1 <- iris[1:10,]
dat2 <- iris[6:15,]

dat1
dat2

nrow(dat1)
nrow(dat2)

# 02 Computing Gower’s distance ----
gower_dist(dat1, dat2)

# If one data frame has less records than the other, the shortest one is 
# recycled over (just like when you’re adding two vectors of unequal length)
gower_dist(iris[1,], dat1)
gower_dist(iris[1:3,], dat1)


# 03 Computing the top-n matches ----
# The function gower_topn returns a list with two arrays.
dat1 <- iris[1:10,]
L <- gower_topn(
  x=dat1, 
  y=iris, 
  n = 3
)


L[[1]]

# The first array is called index. Each column corresponds to one row of x.

# The entries of each column index the top n best matches of that row in x 
# with rows in y

# In this example, the best match of the first row of dat1 is record number 1 
# from iris (this should be obvious, since they are the same record). 
# The second best match is record number 18 from iris.


L[[2]]
# The second array is called distance and it contains the corresponding 
# distances.

L[[2]][1][1]
# It is de distance between observation 1 in x and observation 1 in y

L[[2]][2][1]
# It is de distance between observation 1 in x and observation 18 in y

# 04 Using cluster package ----
data(flower)
head(flower)
nrow(flower)

str(flower)
levels(flower$V4)

df1 = flower[1,]
df2 = flower[2,]

df_total = rbind(
  df1,
  df2
)

?gower_dist

gower_distance = daisy(
  df_total,
  metric = "gower"
)

dd <- as.matrix(
  gower_distance
)

dd


gower_dist(
  df1,
  df2
)

# Both functions return different value for the distance between the first
# and the second observation. This is probably because variables V1 - V6 are
# factors and function gower_dist is treating as continuous variables. Let's
# check it

# 999 Trash ----
if(FALSE){
  ?gower_dist
}
