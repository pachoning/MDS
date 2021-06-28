load_libraries <- function(packages_list){
  
  for(pkg in packages_list){
    print(pkg)
    if(!require(pkg, character.only = TRUE)){
      message(paste0("Installing", pkg, "------------------------------"))
      install.packages(eval(pkg))
      library(pkg, character.only=TRUE)
    }
  }
  
  
}

load_libraries(packages_list=c("pdist", "sets", "stringi", "parallel"))
