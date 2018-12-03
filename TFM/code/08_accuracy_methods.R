# ftp://ftp.auckland.ac.nz/pub/software/CRAN/doc/packages/micEcdat.pdf
library(Ecdat)
source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")
source("tools/fast_mds.R")
source("tools/compute_accuracy.R")
source("tools/build_random_matrix.R")

load("data/Bike-Sharing-Dataset/df_split.RData")


# Accuracy
data(BudgetFood)
x = BudgetFood[1:1000,] %>% select(-sex, -town) 

x = as.data.frame(matrix(rnorm(5*3*10^3), ncol = 5))

x = as.data.frame(
  build_matrix(
    n = 3*10^3,
    p = 5,
    corr_coef = 0
  )
)
dim(x)
# Params for divide and conquer
n_groups = 20

# Params for fast
s = 2
l = 1000
k = 3

metric = "euclidean"
# Classical
if(FALSE){
  starting_time = proc.time()
  rm(results_classical_mds)
  results_classical_mds = classical_mds(
    x = x,
    number_coordinates = s,
    metric = metric
  )
  
  diff_time = proc.time() - starting_time 
  message(paste0("Elapsed time for classical: ", round(diff_time[3], 4), " seconds" ))
}



# Divide and conquer
if(FALSE){
  starting_time = proc.time()
  rm(mds_divide_conquer)
  mds_divide_conquer = divide_conquer_mds(
    x = x,
    groups =  sample(x = n_groups, size = nrow(x), replace = TRUE),
    number_coordinates = s,
    metric = metric
  )
  diff_time = proc.time() - starting_time 
  message(paste0("Elapsed time for divide and conquer: ", round(diff_time[3], 4), " seconds" ))
}


# Fast
if(FALSE){
  starting_time = proc.time()
  rm(mds_fast)
  mds_fast = fast_mds(
    x = x,
    n = nrow(x),
    l = l,
    s = s,
    k = k,
    metric = metric
  )
  
  diff_time = proc.time() - starting_time 
  message(paste0("Elapsed time for fast: ", round(diff_time[3], 4), " seconds" ))
}


# Align divide and classical
results_compare_divide_conquer = compare_methods(
  mds_new_approach = mds_divide_conquer,
  mds_classical = results_classical_mds
)


# Comparing coordinates for divide and classical
plot(mds_divide_conquer[,1], results_compare_divide_conquer$mds_classical_transformed[,1])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_divide_conquer[,2], results_compare_divide_conquer$mds_classical_transformed[,2])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_divide_conquer[,3], results_compare_divide_conquer$mds_classical_transformed[,3])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_divide_conquer[,4], results_compare_divide_conquer$mds_classical_transformed[,4])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_divide_conquer[,5], results_compare_divide_conquer$mds_classical_transformed[,5])
abline(a = 0, b = 1, col = 2, lwd = 2)


# Align fast and classical
results_compare_fast = compare_methods(
  mds_new_approach = mds_fast,
  mds_classical = results_classical_mds
)


# Comparing coordinates for fast and classical
plot(mds_fast[,1], results_compare_fast$mds_classical_transformed[,1])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_fast[,2], results_compare_fast$mds_classical_transformed[,2])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_fast[,3], results_compare_fast$mds_classical_transformed[,3])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_fast[,4], results_compare_fast$mds_classical_transformed[,4])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_fast[,5], results_compare_fast$mds_classical_transformed[,5])
abline(a = 0, b = 1, col = 2, lwd = 2)


# Comparing coordinates for fast and divide
results_compare_divide_conquer_fast = compare_methods(
  mds_new_approach = mds_divide_conquer,
  mds_classical = mds_fast
)

# Comparing coordinates for divide and classical
plot(mds_divide_conquer[,1], results_compare_divide_conquer_fast$mds_classical_transformed[,1])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_divide_conquer[,2], results_compare_divide_conquer_fast$mds_classical_transformed[,2])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_divide_conquer[,3], results_compare_divide_conquer_fast$mds_classical_transformed[,3])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_divide_conquer[,4], results_compare_divide_conquer_fast$mds_classical_transformed[,4])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_divide_conquer[,5], results_compare_divide_conquer_fast$mds_classical_transformed[,5])
abline(a = 0, b = 1, col = 2, lwd = 2)



row.names(mds_fast)
plot(results_classical_mds)
cov(results_classical_mds)

plot(mds_fast[,1], results_classical_mds[,1])
