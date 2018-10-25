library(pracma)


##  Procrustes
U <- randortho(5) # random orthogonal matrix
I = eye(5)
P <- procrustes(U, I)

norm(U-I%*%P$Q, type = 'F')
P$d

