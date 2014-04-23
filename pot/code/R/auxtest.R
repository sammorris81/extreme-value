library(fields)
library(geoR)

s1 <- seq(1, 3)
s2 <- seq(1, 3)
s <- expand.grid(s1, s2)
ns <- nrow(s)

nt <- 3
nknots <- 2
knots <- vector(mode="list", length = nt)
for(t in 1:nt){
	knots[[t]] <- matrix(runif(nknots*2, 1, 3), nknots, 2)
}

membership <- Membership(s, knots)

X <- array(rnorm(nrow(s)*nt*4), dim=c(nrow(s), nt, 4))
X[, , 1] <- 1
beta <- seq(1:4)

  Xbeta <- matrix(NA, ns, nt)
  for (t in 1:nt) {
  	Xbeta[, t] <- X[, t, ] %*% beta
  }
