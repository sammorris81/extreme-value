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
legend1 <- as.expression(bquote(alpha==0.0))
legend2 <- as.expression(bquote(alpha==0.2))
legend3 <- as.expression(bquote(alpha==0.5))
legend4 <- as.expression(bquote(alpha==0.8))
legend5 <- as.expression(bquote(alpha==1.0))
legend("bottomright", legend=c(legend1, legend2, legend3, legend4, legend5), lty=1,
       col=c("black", "firebrick3", "dodgerblue3", "darkolivegreen3", "orange3"))
