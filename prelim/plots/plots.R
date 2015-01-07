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