library(MCMCpack)
A = matrix(c(1,1,-1,0, 1,4), nrow = 3)
p = 2
rotation_matrix = matrix(c(cos(210), sin(210), -sin(210), cos(210)), nrow = 2)
b = c(1,2)
ones_vector = rep(1, nrow(A))

B = ones_vector %*% t(b)

A_modified = p*A%*%rotation_matrix + B

proc = smacof::Procrustes(
  X = A, 
  Y = A_modified
)

proc$rotation
rotation_matrix

1/proc$dilation
p


proc$translation
b 

translation = A - 1/proc$dilation*A_modified%*%proc$rotation


vegan_proc = vegan::procrustes(
  X = A,
  Y = A_modified
)

vegan_proc$rotation
rotation_matrix

vegan_proc$scale
p



vegan_proc$translation
b 

mcm_proc = MCMCpack::procrustes(
  X = A, 
  Xstar = A_modified, 
  translation=TRUE, 
  dilation=TRUE
)

mcm_proc$X.new

ones_vector = rep(1, nrow(A))
mcm_proc$s*A%*%mcm_proc$R + ones_vector %*% t(mcm_proc$tt)
A_modified
