rm(list=ls())
set.seed(20)
source('auxfunctions.R')
source('ts_knots.R')
source('ts_tau.R')
source('ts_z.R')

# need nt, nknots, s, phi.w, tau.alpha, tau.beta, phi.tau, tau, phi.z
set.seed(20)
s <- cbind(runif(15, 0, 10), runif(15, 0, 10))
tau.alpha <- 3
tau.beta  <- 8
phi.w.t   <- 0.9
phi.tau.t <- 0.8
phi.z.t   <- 0.6
nknots  <- 3
nt      <- 10
ns      <- nrow(s)

phi.w.cover <- 0
for (d in 1:100) {
  set.seed(d)
  y <- matrix(0, ns, nt)
  knots <- makeknotsTS(nt, nknots, s, phi.w.t)
  knots.temp <- array(0, dim=c(nknots, 2, nt))
  x.range <- max(s[, 1]) - min(s[, 1])
  y.range <- max(s[, 2]) - min(s[, 2])
  s.range <- max(x.range, y.range)
  knots.temp[, 1, ] <- (knots[, 1, ] - min(s[, 1])) / s.range
  knots.temp[, 2, ] <- (knots[, 2, ] - min(s[, 2])) / s.range
  knots.con <- qnorm(knots.temp)
  phi.w.keep <- rep(NA, 5000)
  phi.w <- att.phi.w <- acc.phi.w <- 0
  mh.phi.w <- 1
  for (iter in 1:5000) {
    knots.mh <- updateKnotsTS(y=y, phi=phi.w, knots=knots, knots.con=knots.con,
                              ts=TRUE, att.phi=att.phi.w, acc.phi=acc.phi.w,
                              mh.phi=mh.phi.w)
    phi.w <- phi.w.keep[iter] <- knots.mh$phi
    acc.phi.w <- knots.mh$acc.phi
    att.phi.w <- knots.mh$att.phi
  }
  # plot(phi.w.keep[10:5000], type="l")
  lower <- quantile(phi.w.keep[10:5000], probs=c(0.025))
  upper <- quantile(phi.w.keep[10:5000], probs=c(0.975))
  # abline(h=c(upper, lower), lty=2)
  if ((phi.w.t < upper) & (phi.w.t > lower)) {
    phi.w.cover <- phi.w.cover + 1
  }
  print(d)
}

phi.tau.cover <- 0
for (d in 1:100) {
  set.seed(d)
  tau <- maketauTS(nt, nknots, tau.alpha, tau.beta, phi.tau.t)
  res <- matrix(0, ns, nt)
  phi.tau.keep <- rep(NA, 5000)
  phi.tau <- att.phi.tau <- acc.phi.tau <- 0
  mh.phi.tau <- 1
  for (iter in 1:5000) {
    tau.mh <- updateTauTS(phi=phi.tau, tau=tau, res=res, ts=TRUE,
                          tau.alpha=tau.alpha, tau.beta=tau.beta,
                          att.phi=att.phi.tau, acc.phi=acc.phi.tau, mh.phi=mh.phi.tau)
    phi.tau <- phi.tau.keep[iter] <- tau.mh$phi
    acc.phi.tau <- tau.mh$acc.phi
    att.phi.tau <- tau.mh$att.phi
  }
  # plot(phi.tau.keep[10:5000], type="l")
  lower <- quantile(phi.tau.keep[10:5000], probs=c(0.025))
  upper <- quantile(phi.tau.keep[10:5000], probs=c(0.975))
  # abline(h=c(upper, lower), lty=2)
  if ((phi.tau.t < upper) & (phi.tau.t > lower)) {
    phi.tau.cover <- 1 + phi.tau.cover
  }
  print(d)
}

phi.z.cover <- 0
for (d in 1:100) {
  set.seed(d)
  z <- makezTS(nt, nknots, tau, phi.z.t, zstar=TRUE)
  tau <- maketauTS(nt, nknots, tau.alpha, tau.beta, phi.tau.t)
  y <- matrix(0, ns, nt)
  phi.z.keep <- rep(NA, 5000)
  phi.z <- att.phi.z <- acc.phi.z <- 0.5
  mh.phi.z <- 1
  for (iter in 1:5000) {
    z.mh <- updateZTS(phi=phi.z, z.star=z, y=y, tau=tau,
                      att.phi=att.phi.z, acc.phi=acc.phi.z, mh.phi=mh.phi.z)
    phi.z <- phi.z.keep[iter] <- z.mh$phi
    acc.phi.z <- tau.mh$acc.phi
    att.phi.z <- tau.mh$att.phi
  }
  # plot(phi.z.keep[10:5000], type="l")
  lower <- quantile(phi.z.keep[10:5000], prob=c(0.025))
  upper <- quantile(phi.z.keep[10:5000], prob=c(0.975))
  # abline(h=c(upper, lower), lty=2)
  if ((phi.z.t < upper) & (phi.z.t > lower)) {
    phi.z.cover <- phi.z.cover + 1
  }
  print(d)
}