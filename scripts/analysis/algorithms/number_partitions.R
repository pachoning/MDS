# Fast
n <- 10^6
s <- 2*10
l <- 800
p_fast <- l/s
n_i <- n/p_fast
i <- 1
while (n_i>l) {
  i <- i + 1
  n_i <- n_i/p_fast
}
p_fast^i
n/p_fast^i
message(paste0("Fast:\nNumber of partitions: ", 
               p_fast^i, 
               "\nSample size of each partition: ", 
               n_i, 
               "\n-------------------------------------------"))

# Divide
p_divide <- 1 + (n-l)/(l-s)
message(paste0("Divide:\nNumber of partitions: ", 
               p_divide, 
               "\nSample size of each partition: ",
               n/p_divide))
