# libraries
library(fields)
library(geoR)

# necessary functions
source('../../../code/R/mcmc.R')
source('../../../code/R/auxfunctions.R')

# data settings
beta.t <- c(10, 2, -3)
nu.t <- 0.5
alpha.t <- 0.9
mixprob.t <- c(0, 1, 1, 1, 1, 0.5)  # 0: Gaussian, 1: t
nknots.t <- c(1, 1, 5, 1, 5, 1)
gau.rho.t <- c(0.10, 0.10, 0.10, 0.10, 0.10, 0.10)
t.rho.t <- c(0.10, 0.10, 0.10, 0.10, 0.10, 0.40)
z.alpha.t <- c(0, 0, 0, 3, 3, 0)
tau.alpha.t <- 2
tau.beta.t  <- 8

# covariate data
s         <- cbind(runif(500), runif(500))
ns        <- nrow(s)
nt        <- 50
nsets     <- 1
nsettings <- length(mixprob.t)
ntest     <- floor(ns / 2)

x <- array(1, c(ns, nt, 3))
for (t in 1:nt) {
    x[, t, 2] <- s[, 1]
    x[, t, 3] <- s[, 2]
}

# Storage for datasets
y <- array(NA, dim=c(ns, nt, nsets, nsettings))
tau.t <- z.t <- knots.t <- vector("list", length=nsettings)

for (setting in 1:nsettings) {
  nknots <- nknots.t[setting]
  tau.t.setting <- array(NA, dim=c(nknots, nt, nsets))
  z.t.setting <- array(NA, dim=c(nknots, nt, nsets))
  knots.t.setting <- array(NA, dim=c(nknots, 2, nt, nsets))
  for (set in 1:nsets) {
    set.seed(setting*100 + set)
    data <- rpotspat(nt=nt, x=x, s=s, beta=beta.t, alpha=alpha.t, nu=nu.t,
                     gau.rho=gau.rho.t[setting], t.rho=t.rho.t[setting],
                     mixprob=mixprob.t[setting], z.alpha=z.alpha.t[setting],
                     tau.alpha=tau.alpha.t, tau.beta=tau.beta.t,
                     nknots=nknots.t[setting])

    y[, , set, setting]        <- data$y
    tau.t.setting[, , set]     <- data$tau
    z.t.setting[, , set]       <- data$z
    knots.t.setting[, , , set] <- data$knots
  }

  tau.t[[setting]]   <- tau.t.setting
  z.t[[setting]]     <- z.t.setting
  knots.t[[setting]] <- knots.t.setting
}

# chi plots
# bin information
d <- as.vector(dist(s))
j <- 1:(nrow(s) - 1)
i <- 2:nrow(s)
ij <- expand.grid(i, j)
ij <- ij[(ij[1] > ij[2]), ]
sites <- cbind(d, ij)
dist <- rdist(s)
diag(dist) <- 0

bins <- seq(0, 1, length=10)
bins <- c(bins, 1.5)

probs <- c(0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)

y.set <- y[, , 1, 1]
threshs <- quantile(y.set, probs=probs, na.rm=T)
exceed.1 <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
for (thresh in 1:length(threshs)){
  exceed.thresh <- att <- acc <- rep(0, (length(bins) - 1))
  for (t in 1:nt) {
    these <- which(y.set[, t] > threshs[thresh])
    n.na <- sum(is.na(y.set[, t]))
    for (b in 1:(length(bins) - 1)) {
      for (site in these){
        inbin <- (dist[site, ] > bins[b]) & (dist[site, ] < bins[b+1])
        att[b] <- att[b] + sum(inbin) - sum(is.na(y.set[inbin, t]))
        acc[b] <- acc[b] + sum(y.set[inbin, t] > threshs[thresh], na.rm=T)
      }
      if (att[b] == 0) {
        exceed.thresh[b] <- 0
      } else {
        exceed.thresh[b] <- acc[b] / att[b]
      }
    }
  }
  exceed.1[, thresh] <- exceed.thresh
}

y.set <- y[, , 1, 2]
threshs <- quantile(y.set, probs=probs, na.rm=T)
exceed.2 <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
for (thresh in 1:length(threshs)){
  exceed.thresh <- att <- acc <- rep(0, (length(bins) - 1))
  for (t in 1:nt) {
    these <- which(y.set[, t] > threshs[thresh])
    n.na <- sum(is.na(y.set[, t]))
    for (b in 1:(length(bins) - 1)) {
      for (site in these){
        inbin <- (dist[site, ] > bins[b]) & (dist[site, ] < bins[b+1])
        att[b] <- att[b] + sum(inbin) - sum(is.na(y.set[inbin, t]))
        acc[b] <- acc[b] + sum(y.set[inbin, t] > threshs[thresh], na.rm=T)
      }
      if (att[b] == 0) {
        exceed.thresh[b] <- 0
      } else {
        exceed.thresh[b] <- acc[b] / att[b]
      }
    }
  }
  exceed.2[, thresh] <- exceed.thresh
}

par(mfrow=c(1, 2))
xplot <- bins[-11]
plot(xplot, exceed[, 1], type="b", ylim=c(0, 0.75), ylab="exceed", xlab="bin distance", pch=1, lty=3, main="chi-plot: t, 1 partition")
for (line in 2:9) {
	lines(xplot, exceed[, line], lty=3)
	points(xplot, exceed[, line], pch=line)
}
legend("topright", lty=1:9, pch=1:9, legend=probs, title="sample quants")


y.set <- y[, , 1, 3]
threshs <- quantile(y.set, probs=probs, na.rm=T)
exceed.3 <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
for (thresh in 1:length(threshs)){
  exceed.thresh <- att <- acc <- rep(0, (length(bins) - 1))
  for (t in 1:nt) {
    these <- which(y.set[, t] > threshs[thresh])
    n.na <- sum(is.na(y.set[, t]))
    for (b in 1:(length(bins) - 1)) {
      for (site in these){
        inbin <- (dist[site, ] > bins[b]) & (dist[site, ] < bins[b+1])
        att[b] <- att[b] + sum(inbin) - sum(is.na(y.set[inbin, t]))
        acc[b] <- acc[b] + sum(y.set[inbin, t] > threshs[thresh], na.rm=T)
      }
      if (att[b] == 0) {
        exceed.thresh[b] <- 0
      } else {
        exceed.thresh[b] <- acc[b] / att[b]
      }
    }
  }
  exceed.3[, thresh] <- exceed.thresh
}

xplot <- bins[-11]
plot(xplot, exceed[, 1], type="b", ylim=c(0, 0.75), ylab="exceed", xlab="bin distance", pch=1, lty=3, main="chi-plot: t, 5 partitions")
for (line in 2:9) {
	lines(xplot, exceed[, line], lty=3)
	points(xplot, exceed[, line], pch=line)
}
legend("topright", lty=1:9, pch=1:9, legend=probs, title="sample quants")


y.set <- y[, , 1, 4]
threshs <- quantile(y.set, probs=probs, na.rm=T)
exceed.4 <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
for (thresh in 1:length(threshs)){
  exceed.thresh <- att <- acc <- rep(0, (length(bins) - 1))
  for (t in 1:nt) {
    these <- which(y.set[, t] > threshs[thresh])
    n.na <- sum(is.na(y.set[, t]))
    for (b in 1:(length(bins) - 1)) {
      for (site in these){
        inbin <- (dist[site, ] > bins[b]) & (dist[site, ] < bins[b+1])
        att[b] <- att[b] + sum(inbin) - sum(is.na(y.set[inbin, t]))
        acc[b] <- acc[b] + sum(y.set[inbin, t] > threshs[thresh], na.rm=T)
      }
      if (att[b] == 0) {
        exceed.thresh[b] <- 0
      } else {
        exceed.thresh[b] <- acc[b] / att[b]
      }
    }
  }
  exceed.4[, thresh] <- exceed.thresh
}

par(mfrow=c(1, 2))
xplot <- bins[-11]
plot(xplot, exceed[, 1], type="b", ylim=c(0, 1), ylab="exceed", xlab="bin distance", pch=1, lty=3, main="chi-plot: skew-t, 1 partition")
for (line in 2:9) {
	lines(xplot, exceed[, line], lty=3)
	points(xplot, exceed[, line], pch=line)
}
legend("topright", lty=1:9, pch=1:9, legend=probs, title="sample quants")


y.set <- y[, , 1, 5]
threshs <- quantile(y.set, probs=probs, na.rm=T)
exceed.5 <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
for (thresh in 1:length(threshs)){
  exceed.thresh <- att <- acc <- rep(0, (length(bins) - 1))
  for (t in 1:nt) {
    these <- which(y.set[, t] > threshs[thresh])
    n.na <- sum(is.na(y.set[, t]))
    for (b in 1:(length(bins) - 1)) {
      for (site in these){
        inbin <- (dist[site, ] > bins[b]) & (dist[site, ] < bins[b+1])
        att[b] <- att[b] + sum(inbin) - sum(is.na(y.set[inbin, t]))
        acc[b] <- acc[b] + sum(y.set[inbin, t] > threshs[thresh], na.rm=T)
      }
      if (att[b] == 0) {
        exceed.thresh[b] <- 0
      } else {
        exceed.thresh[b] <- acc[b] / att[b]
      }
    }
  }
  exceed.5[, thresh] <- exceed.thresh
}

xplot <- bins[-11]
plot(xplot, exceed[, 1], type="b", ylim=c(0, 1), ylab="exceed", xlab="bin distance", pch=1, lty=3, main="chi-plot: skew-t, 5 partitions")
for (line in 2:9) {
	lines(xplot, exceed[, line], lty=3)
	points(xplot, exceed[, line], pch=line)
}
legend("topright", lty=1:9, pch=1:9, legend=probs, title="sample quants")

# lty: non-skew vs skew
# pch: gaussian vs t
# color: 1 partition vs 5 partitions
methods <- c("Gaussian", "t, K=1", "t, K=5", "skew-t, K=1", "skew-t, K=5")
bg <- c("firebrick1", "firebrick1", "dodgerblue1", "firebrick1", "dodgerblue1")
col <- c("firebrick4", "firebrick4", "dodgerblue4", "firebrick4", "dodgerblue4")
pch <- c(24, 22, 22, 22, 22)
lty <- c(1, 1, 1, 3, 3)

par(mfrow=c(2, 2))
# chi-plot sample quantile 0.90
xplot <- bins[-11]
plot(xplot, exceed.1[, 3], type="o", ylim=c(0, 1), ylab="exceed", xlab="bin distance", pch=pch[1], lty=lty[1], col=col[1], bg=bg[1], main="sample quantile 0.90")
lines(xplot, exceed.2[, 3], type="o", lty=lty[2], pch=pch[2], col=col[2], bg=bg[2])
lines(xplot, exceed.3[, 3], type="o", lty=lty[3], pch=pch[3], col=col[3], bg=bg[3])
lines(xplot, exceed.4[, 3], type="o", lty=lty[4], pch=pch[4], col=col[4], bg=bg[4])
lines(xplot, exceed.5[, 3], type="o", lty=lty[5], pch=pch[5], col=col[5], bg=bg[5])
abline(h=0.10, lty=2)
# legend("topright", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg)

# chi-plot sample quantile 0.95
xplot <- bins[-11]
plot(xplot, exceed.1[, 4], type="o", ylim=c(0, 1), ylab="exceed", xlab="bin distance", pch=pch[1], lty=lty[1], col=col[1], bg=bg[1], main="sample quantile 0.95")
lines(xplot, exceed.2[, 4], type="o", lty=lty[2], pch=pch[2], col=col[2], bg=bg[2])
lines(xplot, exceed.3[, 4], type="o", lty=lty[3], pch=pch[3], col=col[3], bg=bg[3])
lines(xplot, exceed.4[, 4], type="o", lty=lty[4], pch=pch[4], col=col[4], bg=bg[4])
lines(xplot, exceed.5[, 4], type="o", lty=lty[5], pch=pch[5], col=col[5], bg=bg[5])
abline(h=0.05, lty=2)
# legend("topright", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg)

# chi-plot sample quantile 0.99
xplot <- bins[-11]
plot(xplot, exceed.1[, 8], type="o", ylim=c(0, 0.7), ylab="exceed", xlab="bin distance", pch=pch[1], lty=lty[1], col=col[1], bg=bg[1], main="sample quantile 0.99")
lines(xplot, exceed.2[, 8], type="o", lty=lty[2], pch=pch[2], col=col[2], bg=bg[2])
lines(xplot, exceed.3[, 8], type="o", lty=lty[3], pch=pch[3], col=col[3], bg=bg[3])
lines(xplot, exceed.4[, 8], type="o", lty=lty[4], pch=pch[4], col=col[4], bg=bg[4])
lines(xplot, exceed.5[, 8], type="o", lty=lty[5], pch=pch[5], col=col[5], bg=bg[5])
abline(h=0.01, lty=2)
# legend("topright", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg)

# chi-plot sample quantile 0.995
xplot <- bins[-11]
plot(xplot, exceed.1[, 9], type="o", ylim=c(0, 0.7), ylab="exceed", xlab="bin distance", pch=pch[1], lty=lty[1], col=col[1], bg=bg[1], main="sample quantile 0.995")
lines(xplot, exceed.2[, 9], type="o", lty=lty[2], pch=pch[2], col=col[2], bg=bg[2])
lines(xplot, exceed.3[, 9], type="o", lty=lty[3], pch=pch[3], col=col[3], bg=bg[3])
lines(xplot, exceed.4[, 9], type="o", lty=lty[4], pch=pch[4], col=col[4], bg=bg[4])
lines(xplot, exceed.5[, 9], type="o", lty=lty[5], pch=pch[5], col=col[5], bg=bg[5])
abline(h=0.005, lty=2)
legend("topright", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg, cex=1.5)

plot(xplot, exceed.1[, 9], type="n", axes=F, xlab="", ylab="")
legend("center", legend=methods, lty=lty, col=col, pch=pch, pt.bg=bg, cex=4, bty="n")

# example of a partition
s1.preds <- seq(1050, 1800, length=30)
s2.preds <- seq(-860, -250, length=30)
s.preds <- expand.grid(s1.preds, s2.preds)
knots <- cbind(runif(5, 1050, 1800), runif(5, -860, -250))
g <- mem(s.preds, knots)
quilt.plot(x=s.preds[, 1], y=s.preds[, 2], z=g, nx=30, ny=30, add.legend=F)
text(knots[1, 1], knots[1, 2], "1", cex=4)
text(knots[2, 1], knots[2, 2], "2", cex=4)
text(knots[3, 1], knots[3, 2], "3", cex=4)
text(knots[4, 1], knots[4, 2], "4", cex=4)
text(knots[5, 1], knots[5, 2], "5", cex=4)
lines(l)

# plot monitoring station ozone locations
load('../../code/analysis/ozone/US-all/us-all-setup.RData')
source('../../code/R/mcmc.R', chdir=T)
source('../../code/R/auxfunctions.R')
plot(s, type="p", xlab="", ylab="")
lines(l)

# plot ozone for days 5 and 34
par(mfrow=c(1, 2))
zlim=range(y[, c(5, 34)], na.rm=T)
quilt.plot(x=s[, 1], y=s[, 2], z=y[, 5], nx=40, ny=40, zlim=zlim, main="Day 5", xlab="", ylab="")
lines(l)
quilt.plot(x=s[, 1], y=s[, 2], z=y[, 34], nx=40, ny=40, zlim=zlim, main="Day 34", xlab="", ylab="")
lines(l)

# chi plot
rm(list=ls())
library(sn)
library(fields)
dest <- function(y, alpha, tau, nu) {
  ft <- dt(y, df=nu, log=T)
  Ft.tau <- pt(tau / sqrt(1 + alpha^2), df=nu, log.p=T)
  Ft <- pt((alpha * y + tau) * sqrt((nu + 1) / (nu + y^2)), df=(nu+1), log.p=T)
  result.log <- ft - Ft.tau + Ft
  return(exp(result.log))
}

d <- seq(0, 5, 0.05)
omega <- exp(-d)

alpha <- 0
nu <- 2
chi.0.2 <- rep(0, length(omega))
chi.0.2[1] <- 1
for (i in 2:length(omega)) {
  omega.i <- omega[i]
  alpha.j <- alpha
  alpha.i <- alpha.j * sqrt(1 - omega.i)
  tau.i   <- sqrt(nu + 1) * (alpha.i + alpha.j * omega.i)
  y.i <- (1 - omega.i) * sqrt(nu + 1) / sqrt(1 - omega.i^2)
  chi.0.2[i] <- 2 * integrate(dest, lower=y.i, upper=Inf,
                          alpha=alpha.i, tau=tau.i, nu=nu)$value
}

alpha <- 10
nu <- 2
chi.10.2 <- rep(0, length(omega))
chi.10.2[1] <- 1
for (i in 2:length(omega)) {
  omega.i <- omega[i]
  alpha.j <- alpha
  alpha.i <- alpha.j * sqrt(1 - omega.i)
  tau.i   <- sqrt(nu + 1) * (alpha.i + alpha.j * omega.i)
  y.i <- (1 - omega.i) * sqrt(nu + 1) / sqrt(1 - omega.i^2)
  chi.10.2[i] <- 2 * integrate(dest, lower=y.i, upper=Inf,
                          alpha=alpha.i, tau=tau.i, nu=nu)$value
}

alpha <- 0
nu <- 3
chi.0.3 <- rep(0, length(omega))
chi.0.3[1] <- 1
for (i in 2:length(omega)) {
  omega.i <- omega[i]
  alpha.j <- alpha
  alpha.i <- alpha.j * sqrt(1 - omega.i)
  tau.i   <- sqrt(nu + 1) * (alpha.i + alpha.j * omega.i)
  y.i <- (1 - omega.i) * sqrt(nu + 1) / sqrt(1 - omega.i^2)
  chi.0.3[i] <- 2 * integrate(dest, lower=y.i, upper=Inf,
                          alpha=alpha.i, tau=tau.i, nu=nu)$value
}

alpha <- 10
nu <- 3
chi.10.3 <- rep(0, length(omega))
chi.10.3[1] <- 1
for (i in 2:length(omega)) {
  omega.i <- omega[i]
  alpha.j <- alpha
  alpha.i <- alpha.j * sqrt(1 - omega.i)
  tau.i   <- sqrt(nu + 1) * (alpha.i + alpha.j * omega.i)
  y.i <- (1 - omega.i) * sqrt(nu + 1) / sqrt(1 - omega.i^2)
  chi.10.3[i] <- 2 * integrate(dest, lower=y.i, upper=Inf,
                          alpha=alpha.i, tau=tau.i, nu=nu)$value
}

alpha <- 0
nu <- 4
chi.0.4 <- rep(0, length(omega))
chi.0.4[1] <- 1
for (i in 2:length(omega)) {
  omega.i <- omega[i]
  alpha.j <- alpha
  alpha.i <- alpha.j * sqrt(1 - omega.i)
  tau.i   <- sqrt(nu + 1) * (alpha.i + alpha.j * omega.i)
  y.i <- (1 - omega.i) * sqrt(nu + 1) / sqrt(1 - omega.i^2)
  chi.0.4[i] <- 2 * integrate(dest, lower=y.i, upper=Inf,
                          alpha=alpha.i, tau=tau.i, nu=nu)$value
}

alpha <- 10
nu <- 4
chi.10.4 <- rep(0, length(omega))
chi.10.4[1] <- 1
for (i in 2:length(omega)) {
  omega.i <- omega[i]
  alpha.j <- alpha
  alpha.i <- alpha.j * sqrt(1 - omega.i)
  tau.i   <- sqrt(nu + 1) * (alpha.i + alpha.j * omega.i)
  y.i <- (1 - omega.i) * sqrt(nu + 1) / sqrt(1 - omega.i^2)
  chi.10.4[i] <- 2 * integrate(dest, lower=y.i, upper=Inf,
                          alpha=alpha.i, tau=tau.i, nu=nu)$value
}

alpha <- 0
nu <- 6
chi.0.6 <- rep(0, length(omega))
chi.0.6[1] <- 1
for (i in 2:length(omega)) {
  omega.i <- omega[i]
  alpha.j <- alpha
  alpha.i <- alpha.j * sqrt(1 - omega.i)
  tau.i   <- sqrt(nu + 1) * (alpha.i + alpha.j * omega.i)
  y.i <- (1 - omega.i) * sqrt(nu + 1) / sqrt(1 - omega.i^2)
  chi.0.6[i] <- 2 * integrate(dest, lower=y.i, upper=Inf,
                          alpha=alpha.i, tau=tau.i, nu=nu)$value
}

alpha <- 10
nu <- 6
chi.10.6 <- rep(0, length(omega))
chi.10.6[1] <- 1
for (i in 2:length(omega)) {
  omega.i <- omega[i]
  alpha.j <- alpha
  alpha.i <- alpha.j * sqrt(1 - omega.i)
  tau.i   <- sqrt(nu + 1) * (alpha.i + alpha.j * omega.i)
  y.i <- (1 - omega.i) * sqrt(nu + 1) / sqrt(1 - omega.i^2)
  chi.10.6[i] <- 2 * integrate(dest, lower=y.i, upper=Inf,
                          alpha=alpha.i, tau=tau.i, nu=nu)$value
}

alpha <- 0
nu <- 20
chi.0.20 <- rep(0, length(omega))
chi.0.20[1] <- 1
for (i in 2:length(omega)) {
  omega.i <- omega[i]
  alpha.j <- alpha
  alpha.i <- alpha.j * sqrt(1 - omega.i)
  tau.i   <- sqrt(nu + 1) * (alpha.i + alpha.j * omega.i)
  y.i <- (1 - omega.i) * sqrt(nu + 1) / sqrt(1 - omega.i^2)
  chi.0.20[i] <- 2 * integrate(dest, lower=y.i, upper=Inf,
                          alpha=alpha.i, tau=tau.i, nu=nu)$value
}

alpha <- 10
nu <- 20
chi.10.20 <- rep(0, length(omega))
chi.10.20[1] <- 1
for (i in 2:length(omega)) {
  omega.i <- omega[i]
  alpha.j <- alpha
  alpha.i <- alpha.j * sqrt(1 - omega.i)
  tau.i   <- sqrt(nu + 1) * (alpha.i + alpha.j * omega.i)
  y.i <- (1 - omega.i) * sqrt(nu + 1) / sqrt(1 - omega.i^2)
  chi.10.20[i] <- 2 * integrate(dest, lower=y.i, upper=Inf,
                          alpha=alpha.i, tau=tau.i, nu=nu)$value
}

chi.gaus <- rep(0, length(omega))

plot(d, chi.0.2, type="l", lty=1, xlim=range(d), ylim=c(0, 1),
  xlab=bquote(h), ylab=bquote(chi(h)),
  main=bquote(paste(chi, " statistic as a function of distance")))
lines(d, chi.10.2, lty=3)
lines(d, chi.0.3, lty=1)
lines(d, chi.10.3, lty=3)
lines(d, chi.0.4, lty=1)
lines(d, chi.10.4, lty=3)
lines(d, chi.0.6, lty=1)
lines(d, chi.10.6, lty=3)
lines(d, chi.0.20, lty=1)
lines(d, chi.10.20, lty=3)
lines(d, chi.gaus, lty=1)

# generate pi(h)
set.seed(200)
rho <- 1
ns <- 400
nsims <- 500
rho <- 1
s <- cbind(runif(ns, 0, 9), runif(ns, 0, 9))
d11 <- rdist(s)

h <- cor <- same.3 <- same.5 <- same.10 <- rep(0, sum(1:(nrow(s) - 1)))
idx <- 0

for (i in 1:(nrow(s) - 1)) {
  cat("Starting i =", i, "\n")
  for (j in ((i+1):nrow(s))) {
    s.ij <- s[c(i, j), ]
    h[idx] <- d11[i, j]
    cor[idx] <- exp(-h[idx])
    for (k in 1:nsims){
      knots <- cbind(runif(3, 0, 9), runif(3, 0, 9))
      d <- rdist(s.ij, knots)
      g <- apply(d, 1, which.min)
      if (g[1] == g[2]) {
        same.3[idx] <- same.3[idx] + 1 / nsims
      }

      knots <- cbind(runif(5, 0, 9), runif(5, 0, 9))
      d <- rdist(s.ij, knots)
      g <- apply(d, 1, which.min)
      if (g[1] == g[2]) {
        same.5[idx] <- same.5[idx] + 1 / nsims
      }

      knots <- cbind(runif(10, 0, 9), runif(10, 0, 9))
      d <- rdist(s.ij, knots)
      g <- apply(d, 1, which.min)
      if (g[1] == g[2]) {
        same.10[idx] <- same.10[idx] + 1 / nsims
      }
    }
    idx <- idx + 1
    if ((j %% 25) == 0) {
      cat("\t j =", j, "\n")
    }
  }
}

d.same <- seq(0, 5, 0.05)
d.window <- 0.025
omega.same <- exp(-d.same)
p.same.3 <- p.same.5 <- p.same.10 <- rep(0, length(omega.same))
p.same.3[1] <- p.same.5[1] <- p.same.10[1] <- 1
for (i in 2:length(omega.same)) {
  d.i <- d.same[i]
  d.l <- d.i - d.window
  d.u <- d.i + d.window
  p.same.3[i] <- mean(same.3[(h >= d.l) & (h < d.u)])
  p.same.5[i] <- mean(same.5[(h >= d.l) & (h < d.u)])
  p.same.10[i] <- mean(same.10[(h >= d.l) & (h < d.u)])
}

save.image(file="pot-chi.RData")

load(file="pot-chi.RData")

plot(d.same, p.same.3, xlim=range(d.same), type="l",
     xlab=bquote(h), ylab=bquote(pi(h)), ylim=c(0, 1),
     main="Simulated probability that two sites will be in the same partition")
lines(d.same, p.same.5, lty=2)
lines(d.same, p.same.10, lty=3)
legend("bottomleft", legend=c("3 knots", "5 knots", "10 knots"), lty=c(1, 2, 3))

# multiply by chi.0.3 and chi10.3 to find partition
chi.0.3.p3 <- p.same.3 * chi.0.3
chi.10.3.p3 <- p.same.3 * chi.10.3
chi.0.3.p5 <- p.same.5 * chi.0.3
chi.10.3.p5 <- p.same.5 * chi.10.3
chi.0.3.p10 <- p.same.10 * chi.0.3
chi.10.3.p10 <- p.same.10 * chi.10.3

# plot with just t-distribution and gaussian
plot(d, chi.0.3, type="l", lty=1, xlim=range(d), ylim=c(0, 1), col="firebrick3",
     xlab=bquote(h), ylab=bquote(chi(h)),
     # main=bquote(paste(chi, " statistic as a function of distance"))
     )
lines(d, chi.10.3, lty=2, col="firebrick3")
lines(d, chi.gaus, lty=1)
legend("topright", lty=c(1, 2, 1),
  legend=c("T", as.expression(bquote(paste("Skew-t, ", alpha==10))), "Gaussian"),
  col=c("firebrick3", "firebrick3", "black"))

plot(d, chi.0.3, type="l", lty=1, xlim=range(d), ylim=c(0, 1), col="firebrick3",
     xlab=bquote(h), ylab=bquote(chi(h)),
     # main=bquote(paste(chi, " statistic as a function of distance"))
     )
lines(d, chi.0.3.p3, lty=1, col="dodgerblue3")
lines(d, chi.0.3.p5, lty=1, col="orange3")
lines(d, chi.0.3.p10, lty=1, col="darkolivegreen3")
lines(d, chi.gaus, lty=1)
legend("topright", lty=c(1, 1, 1, 1, 1),
  legend=c("T, K=1", "T, K=3", "T, K=5", "T, K=10", "Gaussian"),
  col=c("firebrick3", "dodgerblue3", "orange3", "darkolivegreen3", "black"))