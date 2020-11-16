if( require(cluster) == FALSE ){
  message("installing cluster library -----------------------------------")
  install.packages('cluster')
  library(cluster)
}

if( require(fields) == FALSE ){
  message("installing fields library -----------------------------------")
  install.packages('fields')
  library(fields)
}

if( require(ggplot2) == FALSE ){
  message("installing ggplot2 library -----------------------------------")
  install.packages('ggplot2')
  library(ggplot2)
}

if( require(MCMCpack) == FALSE ){
  message("installing MCMCpack library -----------------------------------")
  install.packages('MCMCpack')
  library(MCMCpack)
}

if( require(pdist) == FALSE ){
  message("installing pdist library -----------------------------------")
  install.packages('pdist')
  library(pdist)
}

if( require(rlist) == FALSE ){
  message("installing rlist library -----------------------------------")
  install.packages('rlist')
  library(rlist)
}

if( require(smacof) == FALSE ){
  message("installing smacof library -----------------------------------")
  install.packages('smacof')
  library(smacof)
}

if( require(tidyverse) == FALSE ){
  message("installing tidyverse library -----------------------------------")
  install.packages('tidyverse')
  library(tidyverse)
}
