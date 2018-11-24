install.packages("R.utils")
library(R.utils)


foo <- function() {
  print("Tic")
  for (kk in 1:100) {
    print(kk)
    Sys.sleep(0.1)
  }
  print("Tac")
}



res <- NULL
ss = tryCatch(
  {
    res <- withTimeout(
      {
        foo()
      }, 
      timeout = 2.08
    )
  },

  TimeoutException = function(ex) {
    message("Timeout. Skipping.")
  }
)



res <- withTimeout({
  foo()
  }, 
  timeout = 0.08, 
  onTimeout = "silent"
)

quantile()
?quantile
