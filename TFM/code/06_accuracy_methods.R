# ftp://ftp.auckland.ac.nz/pub/software/CRAN/doc/packages/micEcdat.pdf
library(Ecdat) 
source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")
source("tools/fast_MDS_eigen.R")
source("tools/compute_accuracy.R")
source("tools/build_random_matrix.R")

load("data/Bike-Sharing-Dataset/df_split.RData")


if(FALSE){
  data(BudgetFood)
  x = BudgetFood %>% slice(1:3000) %>% select(-sex, -town) 
  head(x)
}

if(FALSE){
  x = as.data.frame(matrix(rnorm(5*3*10^3), ncol = 5))
}

if(FALSE){
  x = as.data.frame(
    build_matrix(
      n = 3*10^3,
      p = 5,
      corr_coef = 0
    )
  )
}


dim(x)

# Params for fast
s = 2
l = 20
k = 3
metric = "euclidean"
threshold_variance_explained = 0.9

# Classical
if(FALSE){
  rm(mds_classical_sol)
  rm(mds_classical)
  
  starting_time = proc.time()
  
  mds_classical_sol = classical_mds(
    x = x,
    s = s,
    metric = metric
  )
  diff_time = proc.time() - starting_time 
  
  mds_classical = mds_classical_sol$points
  
  
  message(paste0("Elapsed time for classical: ", round(diff_time[3], 4), " seconds" ))
}



# Divide and conquer
if(FALSE){
  2*nrow(x)/l
  rm(mds_divide_conquer_sol)
  rm(mds_divide_conquer)
  
  starting_time = proc.time()
  mds_divide_conquer_sol = divide_conquer_mds(
    x = x,
    l = l,
    s = s,
    metric = metric
  )
  diff_time = proc.time() - starting_time 
  mds_divide_conquer = mds_divide_conquer_sol$points
  
  message(paste0("Elapsed time for divide and conquer: ", round(diff_time[3], 4), " seconds" ))
}



# Fast
if(FALSE){
  
  rm(mds_fast_sol)
  rm(mds_fast)
  
  starting_time = proc.time()
  
  mds_fast_sol = fast_eigen_mds(
    x = x,
    n = nrow(x),
    l = l,
    s = s,
    k = k,
    metric = metric,
    threshold_variance_explained = threshold_variance_explained
  )
  diff_time = proc.time() - starting_time 
  
  mds_fast = mds_fast_sol$points

  message(paste0("Elapsed time for fast: ", round(diff_time[3], 4), " seconds" ))
}


# Align divide and classical
results_compare_divide_conquer = compare_methods(
  mds_new_approach = mds_divide_conquer,
  mds_classical = mds_classical
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
  mds_classical = mds_classical
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

# Comparing coordinates for divide and fast
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



plot(mds_fast[,1], mds_classical[,1])
plot(mds_fast[,2], mds_classical[,2])
plot(mds_fast[,3], mds_classical[,3])
plot(mds_fast[,4], mds_classical[,4])
plot(mds_fast[,5], mds_classical[,5])