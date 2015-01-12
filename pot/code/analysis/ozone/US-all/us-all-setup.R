library(fields)
library(mvtnorm)
library(sn)

rm(list=ls())
source('../../../R/mcmc.R', chdir=T)
source('../../../R/auxfunctions.R')

# Setup from Brian
load("../ozone_data.RData")
S <- cbind(x[s[, 1]], y[s[, 2]])  # expands the grid of x, y locs where we have CMAQ

# Exclude if site is missing more than 50% of its days
excl   <- which(rowMeans(is.na(Y)) > 0.50)
index  <- index[-excl]
Y      <- Y[-excl, ]
S      <- S[-excl, ]

image.plot(x, y, matrix(CMAQ[, 5], nx, ny), main="CMAQ output (ppb) - 2005-07-01 - no thin AQS sites")
points(S)       # Locations of monitoring stations
lines(borders)  # Add state lines

# for presentation
quartz(width=8, height=6)
plot(S, axes=F, ylab="", xlab="", ylim=c(-1600, 1300), xlim=c(-2300, 2400))
lines(borders)

#### start data preprocessing
# only include sites from the eastern US
# keep.these <- which(S[, 1] > 0)
# index      <- index[keep.these]
# Y          <- Y[keep.these, ]
# S          <- S[keep.these, ]
# cmaq       <- CMAQ[index, ]

# take a stratified sample of 800 with equal representation from each region of
# the US
set.seed(4218)
total <- 800
NE <- S[, 1] > 0 & S[, 2] > 0
NW <- S[, 1] < 0 & S[, 2] > 0
SE <- S[, 1] > 0 & S[, 2] < 0
SW <- S[, 1] < 0 & S[, 2] < 0
nNE <- round(mean(NE) * total)
nNW <- round(mean(NW) * total)
nSE <- round(mean(SE) * total)
nSW <- round(mean(SW) * total)
keep.these <- c(sample(which(NE), nNE), sample(which(NW), nNW),
                sample(which(SE), nSE), sample(which(SW), nSW))

index      <- index[keep.these]
Y          <- Y[keep.these, ]
S          <- S[keep.these, ]
cmaq       <- CMAQ[index, ]


# center and scale CMAQ data
cmaq <- (cmaq - mean(cmaq)) / sd(cmaq)

# add in intercept
nt <- ncol(Y)
ns <- nrow(Y)
X <- array(1, dim=c(nrow(cmaq), nt, 2))
for (t in 1:nt) {
  X[, t, 2] <- cmaq[, t]
}

S <- S / 1000
x <- x / 1000
y <- y / 1000

#### end data preprocessing

# Plot CMAQ
image.plot(x, y, matrix(CMAQ[, 5], nx, ny), main="CMAQ output (ppb) - US")
points(S)       # Locations of monitoring stations
lines(borders/1000)  # Add state lines

# Plot a day's data
quilt.plot(x=S[, 1], y=S[, 2], z=Y[, 10], nx=100, ny=100,
           xaxt="n", xlim=c(-2.5, 2.5),
           yaxt="n", ylim=c(-1.65, 1.3),
           main="Ozone values on 10 July 2005")
lines(borders/1000)

# Plot a day's data - presentation
quartz(width=8, height=6)
quilt.plot(x=S[, 1], y=S[, 2], z=Y[, 10], nx=100, ny=100,
           xaxt="n", xlim=c(-2.5, 2.5),
           yaxt="n", ylim=c(-1.65, 1.3),)
lines(borders/1000)

#### 2-fold cross validation
set.seed(1)
cv.idx <- sample(nrow(S), nrow(S), replace=F)
cv.1 <- cv.idx[1:400]
cv.2 <- cv.idx[401:800]
cv.lst <- list(cv.1=cv.1, cv.2=cv.2)

beta.init <- 0
tau.init <- 1

save.image(file="cv-setup-us-all.RData")

# Remove non-US locations (for making maps - not done yet)
east <- s.p[, 1] > 2350
west <- s.p[, 1] < -2250
north <- s.p[, 2] > 1275
south <- s.p[, 2] < -1600

remove <- east | west | north | south


# setting 1: gaussian
# setting 2: t1
# setting 3: t5
# setting 4: t1 T = 0.9
# setting 5: t5 T = 0.9
# setting 6: skew-t1
# setting 7: skew-t1 T = 0.9
# setting 8: skew-t5

# setting 9 - 16: No CMAQ

# chi plot
# bin information
# d <- as.vector(dist(S))
# j <- 1:734
# i <- 2:735
# ij <- expand.grid(i, j)
# ij <- ij[(ij[1] > ij[2]), ]
# sites <- cbind(d, ij)

dist <- rdist(S)
diag(dist) <- 0

bins <- seq(0, 4, 0.5)

probs <- c(0.9, 0.95, 0.99)
threshs <- quantile(Y, probs=probs, na.rm=T)
exceed <- matrix(NA, nrow=(length(bins) - 1), ncol=length(threshs))
for (thresh in 1:length(threshs)){
  exceed.thresh <- att <- acc <- rep(0, (length(bins) - 1))
  for (t in 1:nt) {
    these <- which(Y[, t] > threshs[thresh])
    n.na <- sum(is.na(Y[, t]))
    for (b in 1:(length(bins) - 1)) {
      for (site in these){
        inbin <- (dist[site, ] > bins[b]) & (dist[site, ] < bins[b+1])
        att[b] <- att[b] + sum(inbin) - sum(is.na(Y[inbin, t]))
        acc[b] <- acc[b] + sum(Y[inbin, t] > threshs[thresh], na.rm=T)
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

# par(mfrow=c(1, 2))
xplot <- (1:(length(bins)-1)) - 0.5
plot(xplot, exceed[, 1], type="o", ylim=c(0, 0.35), ylab="exceed", xaxt="n", xlab="bin distance / 1000", pch=1, lty=3, main="chi-plot ozone")
axis(1, at=0:(length(bins)-2), labels=bins[-length(bins)])
for (line in 2:3) {
	lines(xplot, exceed[, line], lty=line)
	points(xplot, exceed[, line], pch=line)
}
legend("topright", lty=1:3, pch=1:3, legend=probs, title="sample quantiles")

# Making the chi plot
# run preprocessing first
# Residuals after lm
y.lm <- Y[, 1]
x.lm <- X[, 1, ]
for (t in 2:nt) {
  y.lm <- c(y.lm, Y[, t])
  x.lm <- rbind(x.lm, X[, t, c(1, 2)])
}
ozone.lm <- lm(y.lm ~ x.lm)
ozone.res <- residuals(ozone.lm)
ozone.int <- ozone.lm$coefficients[1]
ozone.cmaq <- ozone.lm$coefficients[3]
ozone.beta <- c(ozone.int, ozone.cmaq)
res <- matrix(NA, ns, nt)
for (t in 1:nt) {
  res[, t] <- Y[, t] - X[, t, c(1, 2)] %*% ozone.beta
}

dist <- rdist(S)
diag(dist) <- 0

# chi(h) plot
bins.h <- seq(0, 3.5, 0.25)

probs <- c(0.9, 0.95, 0.99)
threshs <- quantile(res, probs=probs, na.rm=T)
exceed.h.res <- matrix(NA, nrow=(length(bins.h) - 1), ncol=length(threshs))

for (thresh in 1:length(threshs)){
  exceed.thresh <- att <- acc <- rep(0, (length(bins.h) - 1))
  for (t in 1:nt) {
    these <- which(res[, t] > threshs[thresh])  # for a given day which sites exceed
    n.na <- sum(is.na(res[, t]))
    for (b in 1:(length(bins.h) - 1)) {
      for (site in these){
        inbin <- (dist[site, ] > bins.h[b]) & (dist[site, ] < bins.h[b+1])
        att[b] <- att[b] + sum(inbin) - sum(is.na(res[inbin, t]))
        acc[b] <- acc[b] + sum(res[inbin, t] > threshs[thresh], na.rm=T)
      }
      if (att[b] == 0) {
        exceed.thresh[b] <- 0
      } else {
        exceed.thresh[b] <- acc[b] / att[b]
      }
    }
  }
  exceed.h.res[, thresh] <- exceed.thresh
}

par(mar=c(5.1, 5.1, 4.1, 2.1))
xplot <- bins.h[-length(bins)] + 0.125
plot(xplot, exceed.h.res[, 1], type="b", pch=21, lty=1, lwd=2,
     xlim=c(0, 3.75), xaxt="n", xlab="bin distance (km) / 1000",
     ylim=c(0, 0.35), ylab=bquote(paste(hat(chi)[c], "(h)")),
     # main=bquote(paste(chi, "(h) for ozone residuals")),
     cex.lab=1.5, cex.axis=1.5, cex.main=2
     )
axis(1, at=bins, cex.axis=1.5)
for (line in 2:3) {
	lines(xplot, exceed.h.res[, line], lty=1, pch=(20 + line), type="b", lwd=2)
}
# abline(v=1, lty=3, lwd=2)
legend("topright", lty=1, pch=21:23, legend=c("c = q(0.90)", "c = q(0.95)", "c = q(0.99)"),
       title="Sample quantiles", cex=1.5, lwd=2, pt.bg="white")

# chi(t) plot
max.lag <- 5
exceed.t.res <- matrix(NA, nrow=max.lag, ncol=length(threshs))

for (thresh in 1:length(threshs)){
  exceed.thresh <- att <- acc <- rep(0, max.lag)
  for (site in 1:ns) {
    these <- which(res[site, ] > threshs[thresh])  # for a site which days exceed
    n.na <- sum(is.na(res[site, ]))
    for (t in these){
      for (lag in 1:max.lag) {
        att[lag] <- att[lag] + 1
        acc[lag] <- acc[lag] + sum(these == (t + lag))
      }
    }
    exceed.thresh <- acc / att
  }
  exceed.t.res[, thresh] <- exceed.thresh
}

xplot <- seq(1, 5, 1)
plot(xplot, exceed.t.res[, 1], type="b", pch=21, lty=1, lwd=2,
     xlim=c(0.5, 5.5), xaxt="n", xlab="lag",
     ylim=c(0, 0.37), ylab=bquote(paste(hat(chi)[c], "(t)")),
     # main=bquote(paste(chi, "(t) for ozone residuals")),
     cex.lab=1.5, cex.axis=1.5, cex.main=2
)
axis(1, at=xplot, cex.axis=1.5)
for (line in 2:3) {
	lines(xplot, exceed.t.res[, line], lty=1, pch=(20 + line), type="b", lwd=2)
}
legend("topright", lty=1, pch=21:23, legend=c("c = q(0.90)", "c = q(0.95)", "c = q(0.99)"),
       title="Sample quantiles", pt.bg="white", cex=1.5, lwd=2)

# qq plot of residuals
library(sn)
res.qq <- sort(as.vector(res))  # sorted residuals
theory.qq <- rep(NA, length(res.qq))
for (i in 1:length(theory.qq)) {
	theory.qq[i] <- i / (length(theory.qq) + 1)
}
res.qq <- (res.qq - mean(res.qq, na.rm=T)) / sd(res.qq)

xplot <- qt(theory.qq, 10)
plot(xplot, res.qq)
abline(0, 1)

xplot <- qst(theory.qq, nu=10, alpha=1)
xplot <- (xplot - mean(xplot))
plot(xplot, res.qq, xlab="Theoretical Quantile", ylab="Observed Quantile", main="Q-Q plot: Skew-t with 10 d.f. and alpha = 1")
abline(0, 1)

# ML - Skew-t fit and diagnostics
library(fields)
library(mvtnorm)
library(sn)

rm(list=ls())
source('../../../R/mcmc.R', chdir=T)
source('../../../R/auxfunctions.R')

# Setup from Brian
load("../ozone_data.RData")
S <- cbind(x[s[, 1]], y[s[, 2]])  # expands the grid of x, y locs where we have CMAQ

# Exclude if site is missing more than 50% of its days
excl  <- which(rowMeans(is.na(Y)) > 0.50)
index <- index[-excl]
Y     <- Y[-excl, ]
S     <- S[-excl, ]
cmaq  <- CMAQ[index, ]
cmaq  <- (cmaq - mean(cmaq)) / sd(cmaq)  # center and scale CMAQ data

# add in intercept
nt <- ncol(Y)
ns <- nrow(Y)
X  <- array(1, dim=c(nrow(cmaq), nt, 2))
for (t in 1:nt) {
  X[, t, 2] <- cmaq[, t]
}

S <- S / 1000
x <- x / 1000
y <- y / 1000

# put all covariate and observations into one long sting
y.lm <- Y[, 1]
x.lm <- X[, 1, ]
for (t in 2:nt) {
  y.lm <- c(y.lm, Y[, t])
  x.lm <- rbind(x.lm, X[, t, c(1, 2)])
}

ozone.stlm  <- selm(y.lm ~ x.lm[, 2], family="ST")  # automatically includes int
ozone.lm    <- lm(y.lm ~ x.lm[, 2])
ozone.res   <- residuals(ozone.stlm)
ozone.dp    <- ozone.stlm@param$dp
ozone.omega <- ozone.dp[3]
ozone.alpha <- ozone.dp[4]
ozone.nu    <- ozone.dp[5]

probs <- seq(0.000001, 0.999999, length=length(ozone.res))
ozone.quant.obs <- quantile(ozone.res, probs=probs)
ozone.quant.the <- rep(0, length(ozone.res))
for (i in 1:33) {
  start <- (i - 1) * 1000 + 1
  if (i < 32) {
    end <- start + 999
  } else {
    end <- length(ozone.res)
  }
  print(paste("start =", start))
  ozone.quant.the[start:end] <- qst(probs[start:end], xi=0, omega=ozone.omega,
                                    alpha=ozone.alpha, nu=ozone.nu)
}
save.image(file="ozone-qq.RData")
plot(ozone.quant.the, ozone.quant.obs, main="Q-Q plot for ozone data")
# set.seed(2087)
# # nw: no obs > 75
# nw.these.l <- which((s[-exceed.75.these, 1] < 0) & (s[-exceed.75.these, 2] >= 0))
# n.nw.these.l <- length(nw.these.l) * 0.25
# nw.thin.l <- sample(nw.these.l, n.nw.these.l, replace=F)

# # nw: at least one obs > 75
# nw.these.h <- which((s[exceed.75.these, 1] < 0) & (s[exceed.75.these, 2] >= 0))
# n.nw.these.h <- length(nw.these.h) * 0.25
# nw.thin.h <- sample(nw.these.h, n.nw.these.h, replace=F)

# # sw: no obs > 75
# sw.these.l <- which((s[-exceed.75.these, 1] < 0) & (s[-exceed.75.these, 2] < 0))
# n.sw.these.l <- length(sw.these.l) * 0.25
# sw.thin.l <- sample(sw.these.l, n.sw.these.l, replace=F)

# # sw: at least one obs > 75
# sw.these.h <- which((s[exceed.75.these, 1] < 0) & (s[exceed.75.these, 2] < 0))
# n.sw.these.h <- length(sw.these.h) * 0.25
# sw.thin.h <- sample(sw.these.h, n.sw.these.h, replace=F)

# # ne: no obs > 75
# ne.these.l <- which((s[-exceed.75.these, 1] >= 0) & (s[-exceed.75.these, 2] >= 0))
# n.ne.these.l <- length(ne.these.l) * 0.25
# ne.thin.l <- sample(ne.these.l, n.ne.these.l, replace=F)

# # ne: at least one obs > 75
# ne.these.h <- which((s[exceed.75.these, 1] >= 0) & (s[exceed.75.these, 2] >= 0))
# n.ne.these.h <- length(ne.these.h) * 0.25
# ne.thin.h <- sample(ne.these.h, n.ne.these.h, replace=F)

# # se: no obs > 75
# se.these.l <- which((s[-exceed.75.these, 1] >= 0) & (s[-exceed.75.these, 2] < 0))
# n.se.these.l <- length(se.these.l) * 0.25
# se.thin.l <- sample(se.these.l, n.se.these.l, replace=F)

# # se: at least one obs > 75
# se.these.h <- which((s[exceed.75.these, 1] >= 0) & (s[exceed.75.these, 2] < 0))
# n.se.these.h <- length(se.these.h) * 0.25
# se.thin.h <- sample(se.these.h, n.se.these.h, replace=F)

# thin     <- c(nw.thin.l, nw.thin.h, sw.thin.l, sw.thin.h, ne.thin.l, ne.thin.h, se.thin.l, se.thin.h)

# index    <- index[thin]
# s        <- cmaq.s[index, ]  # only include the stratified sample of sites
# aqs.cmaq <- CMAQ[thin, ]  # CMAQ measurements at the AQS sites
# aqs      <- Y[thin, ]     # AQS measurements at the AQS sites


# Covariates including lat, long, and CMAQ
# ns <- nrow(aqs)
# nt <- ncol(aqs)
# X <- array(1, dim=c(ns, nt, 7))
# for(t in 1:nt){
  # X[, t, 2] <- s[, 1]    # Long
  # X[, t, 3] <- s[, 2]    # Lat
  # X[, t, 4] <- s[, 1]^2  # Long^2
  # X[, t, 5] <- s[, 2]^2  # Lat^2
  # X[, t, 6] <- s[, 1] * s[, 2]  # Interaction
  # X[, t, 7] <- aqs.cmaq[, t]   # CMAQ
# }

# X.scale <- array(1, dim=c(ns, nt, 7))
# for(t in 1:nt){
  # X[, t, 2] <- s.scale[, 1]    # Long
  # X[, t, 3] <- s.scale[, 2]    # Lat
  # X[, t, 4] <- s.scale[, 1]^2  # Long^2
  # X[, t, 5] <- s.scale[, 2]^2  # Lat^2
  # X[, t, 6] <- s.scale[, 1] * s.scale[, 2]  # Interaction
  # X[, t, 7] <- aqs.cmaq[, t]   # CMAQ
# }


#### 5-fold cross validation
# cv.idx <- sample(nrow(s.scale), nrow(s.scale), replace=F)

# cv.1 <- cv.idx[1:108]
# cv.2 <- cv.idx[108:216]
# cv.3 <- cv.idx[217:324]
# cv.4 <- cv.idx[325:432]
# cv.5 <- cv.idx[433:544]

# cv.lst <- list(cv.1=cv.1, cv.2=cv.2, cv.3=cv.3, cv.4=cv.4, cv.5=cv.5)

# Rescale s so each dimension is in (0, 1)
# S.scale      <- matrix(NA, nrow=nrow(S), ncol=ncol(S))
# S.scale[, 1] <- (S[, 1] - range(S[, 1])[1]) / (range(S[, 1])[2] - range(S[, 1])[1])
# S.scale[, 2] <- (S[, 2] - range(S[, 2])[1]) / (range(S[, 2])[2] - range(S[, 2])[1])

# # Statified sample of 50% to make covariance calculations easier
# # stratified by whether exceed 75
# # which sites have at no days that exceed 75, gives the row number of the index
# gthan.75.these <- which(rowMeans(aqs >= 75, na.rm=T) > 0)
# lthan.75.these <- which(rowMeans(aqs >= 75, na.rm=T) ==0)
# n.gthan.75 <- round(length(gthan.75.these) * 0.5)
# n.lthan.75 <- round(length(lthan.75.these) * 0.5)

# set.seed(2087)
# keep.gthan.these <- sample(gthan.75.these, n.gthan.75, replace=F)
# keep.lthan.these <- sample(lthan.75.these, n.lthan.75, replace=F)
# keep.these <- c(keep.gthan.these, keep.lthan.these)
