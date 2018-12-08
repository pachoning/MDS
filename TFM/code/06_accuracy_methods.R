# ftp://ftp.auckland.ac.nz/pub/software/CRAN/doc/packages/micEcdat.pdf
library(Ecdat) 
source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")
source("tools/fast_MDS_eigen.R")
source("tools/classic_mds.R")
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

# Classic
if(FALSE){
  rm(mds_classic_sol)
  rm(mds_classic)
  
  starting_time = proc.time()
  
  mds_classic_sol = classic.mds(
    x = x,
    s = s,
    metric = metric
  )
  diff_time = proc.time() - starting_time 
  
  mds_classic = mds_classic_sol$points
  
  
  message(paste0("Elapsed time for classic: ", round(diff_time[3], 4), " seconds" ))
}



# Divide and conquer
if(FALSE){
  2*nrow(x)/l
  rm(mds_divide_conquer_sol)
  rm(mds_divide_conquer)
  
  starting_time = proc.time()
  mds_divide_conquer_sol = divide_conquer.mds(
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
  
  mds_fast_sol = fast.mds(
    x = x,
    n = nrow(x),
    l = l,
    s = s,
    k = k,
    metric = metric
  )
  diff_time = proc.time() - starting_time 
  
  mds_fast = mds_fast_sol$points

  message(paste0("Elapsed time for fast: ", round(diff_time[3], 4), " seconds" ))
}


# Align divide and classic
results_compare_divide_conquer = compare.methods(
  mds_new_approach = mds_divide_conquer,
  mds_classic = mds_classic
)


# Comparing coordinates for divide and classic
plot(mds_divide_conquer[,1], results_compare_divide_conquer$mds_classic_transformed[,1])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_divide_conquer[,2], results_compare_divide_conquer$mds_classic_transformed[,2])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_divide_conquer[,3], results_compare_divide_conquer$mds_classic_transformed[,3])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_divide_conquer[,4], results_compare_divide_conquer$mds_classic_transformed[,4])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_divide_conquer[,5], results_compare_divide_conquer$mds_classic_transformed[,5])
abline(a = 0, b = 1, col = 2, lwd = 2)


# Align fast and classic
results_compare_fast = compare.methods(
  mds_new_approach = mds_fast,
  mds_classic = mds_classic
)


# Comparing coordinates for fast and classic
plot(mds_fast[,1], results_compare_fast$mds_classic_transformed[,1])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_fast[,2], results_compare_fast$mds_classic_transformed[,2])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_fast[,3], results_compare_fast$mds_classic_transformed[,3])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_fast[,4], results_compare_fast$mds_classic_transformed[,4])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_fast[,5], results_compare_fast$mds_classic_transformed[,5])
abline(a = 0, b = 1, col = 2, lwd = 2)


# Comparing coordinates for fast and divide
results_compare_divide_conquer_fast = compare.methods(
  mds_new_approach = mds_divide_conquer,
  mds_classic = mds_fast
)

# Comparing coordinates for divide and fast
plot(mds_divide_conquer[,1], results_compare_divide_conquer_fast$mds_classic_transformed[,1])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_divide_conquer[,2], results_compare_divide_conquer_fast$mds_classic_transformed[,2])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_divide_conquer[,3], results_compare_divide_conquer_fast$mds_classic_transformed[,3])
abline(a = 0, b = 1, col = 2, lwd = 2)

plot(mds_divide_conquer[,4], results_compare_divide_conquer_fast$mds_classic_transformed[,4])
abline(a = 0, b = 1, col = 2, lwd = 2)


plot(mds_divide_conquer[,5], results_compare_divide_conquer_fast$mds_classic_transformed[,5])
abline(a = 0, b = 1, col = 2, lwd = 2)



plot(mds_fast[,1], mds_classic[,1])
plot(mds_fast[,2], mds_classic[,2])
plot(mds_fast[,3], mds_classic[,3])
plot(mds_fast[,4], mds_classic[,4])
plot(mds_fast[,5], mds_classic[,5])

