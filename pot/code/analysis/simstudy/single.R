# test single

tic <- proc.time()
for(i in 1:100000){
  X <- matrix(rnorm(10000), 100, 100)
  xtx <- X %*% X
}
toc <- proc.time()

elapsed <- toc[3] - tic[3]
cat("this took", elapsed, "seconds")
