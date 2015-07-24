rm(list=ls())
library(fields)
s <- matrix(c(-1, 0, 1, 0), 2, 2, byrow = TRUE)
rdist(s)

nsamp <- 20000
knots <- array(NA, dim = c(2, 2, nsamp))
same <- rep(NA, nsamp)

set.seed(1)
for (i in 1:nsamp) {
  theta.1 <- runif(1, 0, 2 * pi)
  theta.2 <- runif(1, 0, 2 * pi)
  r.1 <- runif(1, 0, 1)
  r.2 <- runif(1, 0, 1)
  knots[1, 1, i] <- r.1 * cos(theta.1)
  knots[1, 2, i] <- r.1 * sin(theta.1)
  knots[2, 1, i] <- r.2 * cos(theta.2)
  knots[2, 2, i] <- r.2 * sin(theta.2)
  d <- rdist(s, knots[, , i])
  same[i] <- which(d[1, ] == min(d[1, ])) == which(d[2, ] == min(d[2, ]))
}

knots <- array(NA, dim = c(3, 2, nsamp))
for (i in 1:nsamp) {
  theta.1 <- runif(1, 0, 2 * pi)
  theta.2 <- runif(1, 0, 2 * pi)
  theta.3 <- runif(1, 0, 2 * pi)
  r.1 <- runif(1, 0, 1)
  r.2 <- runif(1, 0, 1)
  r.3 <- runif(1, 0, 1)
  knots[1, 1, i] <- r.1 * cos(theta.1)
  knots[1, 2, i] <- r.1 * sin(theta.1)
  knots[2, 1, i] <- r.2 * cos(theta.2)
  knots[2, 2, i] <- r.2 * sin(theta.2)
  knots[3, 1, i] <- r.3 * cos(theta.3)
  knots[3, 2, i] <- r.3 * sin(theta.3)
  d <- rdist(s, knots[, , i])
  same[i] <- which(d[1, ] == min(d[1, ])) == which(d[2, ] == min(d[2, ]))
}


plot(0, 0, ylim=c(-1, 1), xlim=c(-1, 1), type = "n")
col <- 1
for (i in 1:15) {
  if (partition.1[i] == partition.2[i]) {
    points(knots[, , i], col = col)
    col <- col + 1
  }
}

plot(0, 0, ylim=c(-1, 1), xlim=c(-1, 1), type = "n")
r <- seq(0.02, 1, by = 0.02)
theta <- seq(0, 2 * pi, by = 0.01)
knots <- matrix(NA, 2, 2)
r.1 <- runif(1, 0, 1)
theta.1 <- runif(1, 0, 2 * pi)
x.o <- r.1 * cos(theta.1)
y.o <- r.1 * sin(theta.1)
knots[1, ] <- c(x.o, y.o)
knots[2, ] <- c(0, 0)  # initial one
d <- rdist(s, knots)
if (which(d[1, ] == min(d[1, ])) == which(d[2, ] == min(d[2, ]))) {
  points(0, 0, pch = 16, col = "dodgerblue4", cex = 0.75)
} else {
  points(0, 0, pch = 16, col = "firebrick4", cex = 0.75)
}
for (i in 1:length(r)) {
  for (j in 1:length(theta)) {
    x <- r[i] * cos(theta[j])
    y <- r[i] * sin(theta[j])
    knots[2, ] <- c(x, y)
    d <- rdist(s, knots)
    if (which(d[1, ] == min(d[1, ])) == which(d[2, ] == min(d[2, ]))) {
      points(x, y, pch = 16, col = "dodgerblue4", cex=0.75)
    } else {
      points(x, y, pch = 16, col = "firebrick4", cex=0.75)
    }
  }
}
points(x.o, y.o, pch = 16, col = "black")
