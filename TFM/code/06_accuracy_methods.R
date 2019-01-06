# ftp://ftp.auckland.ac.nz/pub/software/CRAN/doc/packages/micEcdat.pdf
source("tools/load_libraries.R")
source("tools/divide_conquer_mds.R")
source("tools/fast_MDS_eigen.R")
source("tools/classical_mds.R")
source("tools/gower_interpolation_mds.R")
source("tools/compute_accuracy.R")
source("tools/build_random_matrix.R")

x = matrix(rnorm(3*1*10^3), ncol = 3)
row.names(x) = 1:nrow(x)
dim(x)

# Params for fast
s = 3
l = 500
k = 3
metric = "euclidean"

# Classical
if(FALSE){
  rm(mds_classical_sol)
  rm(mds_classical)
  
  starting_time = proc.time()
  
  mds_classical_sol = classical.mds(
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
    l = l,
    s = s,
    k = k,
    metric = metric
  )
  diff_time = proc.time() - starting_time 
  
  mds_fast = mds_fast_sol$points

  message(paste0("Elapsed time for fast: ", round(diff_time[3], 4), " seconds" ))
}

# Gower interpolation
if(FALSE){
  
  rm(mds_gower_sol)
  rm(mds_gower)
  
  starting_time = proc.time()
  
  mds_gower_sol = gower.interpolation.mds(
    x = x,
    l = l,
    s = s
  )
  diff_time = proc.time() - starting_time 
  
  mds_gower = mds_gower_sol$points
  
  message(paste0("Elapsed time for fast: ", round(diff_time[3], 4), " seconds" ))
}

# Align divide and classical
results_compare_divide_conquer = compare.methods(
  mds_new_approach = mds_divide_conquer,
  mds_classical = x
)

# Comparing coordinates for divide and classical
pdf(
  file.path(
    getwd(),
    "thesis",
    "images",
    "first_div.pdf"
  ),
  width = 4, 
  height = 4
)
ggplot(
  data.frame(
    x = mds_divide_conquer[,1],
    y = results_compare_divide_conquer$mds_classical_transformed[,1]
  ), 
  aes(x=x, y=y)
) +
  geom_point(shape = 1)  +
  geom_abline(aes(intercept = 0, slope = 1, colour = "B")) + 
  guides(fill=FALSE, color=FALSE) + 
  xlim(-4, 4) + 
  ylim(-4,4)+
  labs(
    x ="First dim. of X", 
    y = "First dim. of Divide and Conquer MDS"
  )
dev.off()

xtable::xtable(
  cor(results_compare_divide_conquer$mds_classical_transformed, mds_divide_conquer)
)


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
results_compare_fast = compare.methods(
  mds_new_approach = mds_fast,
  mds_classical = x
)


# Comparing coordinates for fast and classical
pdf(
  file.path(
    getwd(),
    "thesis",
    "images",
    "third_fast.pdf"
  ),
  width = 4, 
  height = 4
)
ggplot(
  data.frame(
    x = mds_fast[,3],
    y = results_compare_fast$mds_classical_transformed[,3]
  ), 
  aes(x=x, y=y)
) +
  geom_point(shape = 1)  +
  geom_abline(aes(intercept = 0, slope = 1, colour = "B")) + 
  guides(fill=FALSE, color=FALSE) + 
  xlim(-4, 4) + 
  ylim(-4,4) +
  labs(
    x ="Third dim. of X", 
    y = "Third dim. of Fast MDS"
  )
dev.off()


xtable::xtable(
  cor(results_compare_fast$mds_classical_transformed, mds_fast)
)


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


# Align gower and classical
results_compare_gower = compare.methods(
  mds_new_approach = mds_gower,
  mds_classical = x
)

# Comparing coordinates for divide and classical
pdf(
  file.path(
    getwd(),
    "thesis",
    "images",
    "third_gower.pdf"
  ),
  width = 4, 
  height = 4
)
ggplot(
  data.frame(
    x = mds_gower[,3],
    y = results_compare_gower$mds_classical_transformed[,3]
  ), 
  aes(x=x, y=y)
) +
  geom_point(shape = 1)  +
  geom_abline(aes(intercept = 0, slope = 1, colour = "B")) + 
  guides(fill=FALSE, color=FALSE) + 
  xlim(-4, 4) + 
  ylim(-4,4) +
  labs(
    x ="Third dim. of X", 
    y = "Third dim. of MDS based on Gower interp."
  )
dev.off()

xtable::xtable(
  cor(results_compare_gower$mds_classical_transformed, mds_gower)
)




# Comparing coordinates for fast and divide
results_compare_divide_conquer_fast = compare.methods(
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

