# libraries
library(fields)
library(geoR)

# necessary functions
source('../../R/mcmc.R')
source('../../R/auxfunctions.R')

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

# chi plot
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

y.set <- y[, , 1, 2]

probs <- c(0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
threshs <- quantile(y.set, probs=probs, na.rm=T)
exceed <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
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
  exceed[, thresh] <- exceed.thresh
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

probs <- c(0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
threshs <- quantile(y.set, probs=probs, na.rm=T)
exceed <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
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
  exceed[, thresh] <- exceed.thresh
}

xplot <- bins[-11]
plot(xplot, exceed[, 1], type="b", ylim=c(0, 0.75), ylab="exceed", xlab="bin distance", pch=1, lty=3, main="chi-plot: t, 5 partitions")
for (line in 2:9) { 
	lines(xplot, exceed[, line], lty=3)
	points(xplot, exceed[, line], pch=line)
}
legend("topright", lty=1:9, pch=1:9, legend=probs, title="sample quants")

# chi plot
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

y.set <- y[, , 1, 4]

probs <- c(0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
threshs <- quantile(y.set, probs=probs, na.rm=T)
exceed <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
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
  exceed[, thresh] <- exceed.thresh
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

probs <- c(0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
threshs <- quantile(y.set, probs=probs, na.rm=T)
exceed <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
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
  exceed[, thresh] <- exceed.thresh
}

xplot <- bins[-11]
plot(xplot, exceed[, 1], type="b", ylim=c(0, 1), ylab="exceed", xlab="bin distance", pch=1, lty=3, main="chi-plot: skew-t, 5 partitions")
for (line in 2:9) { 
	lines(xplot, exceed[, line], lty=3)
	points(xplot, exceed[, line], pch=line)
}
legend("topright", lty=1:9, pch=1:9, legend=probs, title="sample quants")

