library(future)
plan(multiprocess)
library(promises)
library(tibble)

a %<-% {
  cat("Future 'a' ...")
  Sys.sleep(2)
  cat("done\n")
  log("b")
}

cat("Waiting for 'a' to be resolved ...\n")
f <- futureOf(a)
count <- 1
resolved(f)
while (!resolved(f)) {
  cat(count, "\n")
  Sys.sleep(0.2)
  count <- count + 1
}
cat("Waiting for 'a' to be resolved ... DONE\n")
a

tryCatch(
  log("b"),
  error = "blabla"
)

?tryCatch

i_row = 10000
x = data.frame(
  x1 = rnorm(i_row),
  x2 = rnorm(i_row),
  x3 = rnorm(i_row),
  x4 = rnorm(i_row),
  x5 = rnorm(i_row),
  x6 = rnorm(i_row)
)

pp %<-% {
  daisy(
    x,
    metric = "euclidean"
  )
}


f <- futureOf(pp)
count <- 1
resolved(f)
s = as.matrix(pp)
runif(1000000)

