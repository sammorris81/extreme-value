rm(list=ls())
set.seed(200)
months <- seq(1:36)
month1 <- rnorm(12, 0.5, 0.5) - (months[1:12] - 6.5)^2 / 50
month2 <- rnorm(12, 1, 0.3) - (months[1:12] - 6.5)^2 / 50
month3 <- rnorm(12, 0.5, 0.5) - (months[1:12] - 6.5)^2 / 50

plot(months, c(month1, month2, month3), ylab="Observed value", xlab="Year",
     axes=F, xlim=c(0, 42))
axis(1, at=c(6.5, 18.5, 30.5), labels=c(1, 2, 3), lwd=0)
points(months[c(10, 18, 27)], c(month1[10], month2[6], month3[3]), pch=16,
       col="firebrick3")
abline(v=c(0, 12.5, 24.5, 37), lty=3)
abline(h=0.63)
text(x=37.5, y=0.65, adj=c(0, 0), labels="Threshold")

set.seed(200)
knots <- runif(10, 1, 9)
xplot <- seq(0.75, 9.25, 0.01)
rho <- 0.5
weight <- matrix(0, nrow=length(knots), ncol=length(xplot))
for (i in 1:length(knots)) {
  weight[i, ] <- exp(-0.5 * (xplot - knots[i])^2 / rho)
}
sumweight <- apply(weight, 2, sum)
for(j in 1:length(xplot)) {
  weight[, j] <- weight[, j] / sumweight[j]
}

intensity <- rgamma(length(knots), 10, 3)
intensity.weight <- matrix(0, nrow=length(knots), ncol=length(xplot))
for (k in 1:length(knots)) {
  intensity.weight[k, ] <- intensity[k] * weight[k, ]
}


value.0.0 <- value.0.2 <- value.0.5 <- value.0.8 <- value.1 <- rep(0, length(xplot))

alpha <- 0.5
for (i in 1:length(xplot)) {
  value.0.5[i] <- (sum((intensity * weight[, i])^(1 / alpha)))^alpha
}

alpha <- 0.8
for (i in 1:length(xplot)) {
  value.0.8[i] <- (sum((intensity * weight[, i])^(1 / alpha)))^alpha
}

alpha <- 0.2
for (i in 1:length(xplot)) {
  value.0.2[i] <- (sum((intensity * weight[, i])^(1 / alpha)))^alpha
}

alpha <- 0
for (i in 1:length(xplot)) {
  value.0.0[i] <- max(intensity * weight[, i])
}

alpha <- 1
for (i in 1:length(xplot)) {
  value.1[i] <-sum((intensity * weight[, i]))
}

quartz(width=12, height=6)
par(mfrow=c(1, 2))
range <- range(intensity.weight, intensity)
plot(xplot, intensity[1] * weight[1, ], type="l", ylim=range,
    ylab="Intensity * Weight", xlab="Location", xaxt="n")
axis(1, at=c(1:9))
for (k in 1:length(knots)) {
  lines(xplot, intensity[k] * weight[k, ])
}
points(knots, intensity, pch=16, col="red")

plot(xplot, value.0.2, type="l", ylim=range,
     xaxt="n", xlab="Location", ylab="Observed value", col="firebrick3")
lines(xplot, value.0.5, col="dodgerblue3")
lines(xplot, value.0.8, col="darkolivegreen3")
lines(xplot, value.0.0, col="black")
lines(xplot, value.1, col="orange3")
axis(1, at=c(1:9))
legend5 <- as.expression(bquote(alpha==0.0))
legend4 <- as.expression(bquote(alpha==0.2))
legend3 <- as.expression(bquote(alpha==0.5))
legend2 <- as.expression(bquote(alpha==0.8))
legend1 <- as.expression(bquote(alpha==1.0))
legend("bottomright", legend=c(legend1, legend2, legend3, legend4, legend5), lty=1,
       col=c("orange3", "darkolivegreen3", "dodgerblue3", "firebrick3", "black"))


# link functions
rm(list=ls())
library(evd)
set.seed(200)
x <- rnorm(70, -1, 2)
xplot <- seq(-10, 10, 0.01)
logit <- exp(xplot) / (1 + exp(xplot))
probit <- pnorm(xplot)
cloglog <- 1 - exp(-exp(xplot))
xi <- 0.4
gev <- 1 - pgev(-xplot, shape=xi)
plot(xplot, logit, type="l")
lines(xplot, probit)
lines(xplot, cloglog)
lines(xplot, gev, lty=2)

exceed <- 2
y.norm <- x + rnorm(length(x), 0, 1)
y.bin.norm <- (y.norm > exceed)

y.gev <- x + rgev(length(x), 0, 1, xi)
y.bin.gev <- (y.gev > exceed)

y.gev.2 <- x + rgev(length(x), 0, 1, -xi)
y.bin.gev.2 <- (y.gev.2 > exceed)

probit.1 <- glm(y.bin.gev ~ x, family=binomial(link="probit"))
probit.1.int <- probit.1$coef[1]
probit.1.slope <- probit.1$coef[2]
probit.1.xbeta <- probit.1.int + probit.1.slope * xplot
logit.1 <- glm(y.bin.gev ~ x, family=binomial)
logit.1.int <- logit.1$coef[1]
logit.1.slope <- logit.1$coef[2]
logit.1.xbeta <- logit.1.int + logit.1.slope * xplot

probit.2 <- glm(y.bin.gev.2 ~ x, family=binomial(link="probit"))
probit.2.int <- probit.2$coef[1]
probit.2.slope <- probit.2$coef[2]
probit.2.xbeta <- probit.2.int + probit.2.slope * xplot
logit.2 <- glm(y.bin.gev.2 ~ x, family=binomial)
logit.2.int <- logit.2$coef[1]
logit.2.slope <- logit.2$coef[2]
logit.2.xbeta <- logit.2.int + logit.2.slope * xplot

yplot.prob.1 <- pnorm(probit.1.xbeta)
yplot.logit.1 <- exp(logit.1.xbeta) / (1 + exp(logit.1.xbeta))
yplot.prob.2 <- pnorm(probit.2.xbeta)
yplot.logit.2 <- exp(logit.2.xbeta) / (1 + exp(logit.2.xbeta))
yplot.gev <- 1 - pgev(exceed - xplot, 0, 1, xi)
yplot.cloglog <- 1 - exp(-exp(-(exceed - xplot)))
yplot.gev.2 <- 1 - pgev(exceed - xplot, 0, 1, -xi)

quartz(width=11, height=5.5)
par(mfrow=c(1, 2))
plot(x, y.bin.gev, ylab="P(Y = 1)", xlab="x", xlim=c(-10, 10))
lines(xplot, yplot.prob.1, lty=2)
lines(xplot, yplot.logit.1, lty=3)
lines(xplot, yplot.gev, lty=1, col="dodgerblue3")
lines(xplot, yplot.gev.2, lty=1, col="firebrick3")
legend3 <- as.expression(bquote(paste("GEV, ", xi==0.4)))
legend4 <- as.expression(bquote(paste("GEV, ", xi==-0.4)))
legend("left", legend=c("Probit", "Logistic", legend3, legend4),
       lty=c(2, 3, 1, 1), col=c("black", "black", "dodgerblue3", "firebrick3"))

plot(x, y.bin.gev.2, ylab="P(Y = 1)", xlab="x", xlim=c(-10, 10))
lines(xplot, yplot.prob.2, lty=2)
lines(xplot, yplot.logit.2, lty=3)
lines(xplot, yplot.gev, lty=1, col="dodgerblue3")
lines(xplot, yplot.gev.2, lty=1, col="firebrick3")
legend3 <- as.expression(bquote(paste("GEV, ", xi==0.4)))
legend4 <- as.expression(bquote(paste("GEV, ", xi==-0.4)))
legend("left", legend=c("Probit", "Logistic", legend3, legend4),
       lty=c(2, 3, 1, 1), col=c("black", "black", "dodgerblue3", "firebrick3"))