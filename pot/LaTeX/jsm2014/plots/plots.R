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
load('../../code/analysis/ozone/SE/cv-setup-se.RData')
source('../../code/R/mcmc.R')
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
