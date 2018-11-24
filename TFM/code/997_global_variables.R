library(rlist)

list_mds <<- list()
list_positions <<- list()
n_recursive_calls <<- 0




split_vector = function(
  x,
  times_called = 0
){
  length_vector = length(x)
  n_divisions = 2
  if(length_vector <= 10){
    message(x)
    list_to_return <<- list.append(list_to_return, x)
  }else{
    groups_division = sample(1:n_divisions, length(x), replace=T)
    unique_groups = unique(groups_division)
    for(k_group in unique_groups){
      observations_group = which(groups_division == k_group)
      x_split = x[observations_group]
      times_called = times_called + 1
      split_vector(
        x_split,
        times_called
      )
    }
  }
}

nx = 1:11
split_vector(nx)
ss = split_vector(nx)
dd = list() 
dd[[5]] = 7
message(dd)
