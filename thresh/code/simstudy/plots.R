#### Quantile scores

library(evd)
rm(list=ls())

set.seed(2087)
xi <- 0.25
sig <- 1
n <- 16
x <- round(rgpd(n, 0, scale=sig, shape=xi), 3)
x <- sort(x)

tau <- seq(0, 1, 0.001)

quant <- sig / xi * ((1-tau)^(-xi) - 1)

par(mfrow=c(4, 4))
par(oma=c(2, 2, 2, 0))
par(mar=c(3, 3, 2, 1))
for (i in 1:n) {
	quant.score <- (quant - x[i]) * ((x[i] < quant) - tau)
	
	plot(tau, quant.score, type="l", main=paste("x = ", x[i]), xlab="", ylab="")

}
mtext("tau", side="1", outer=T)
mtext("quantile scores", side="2", outer=T)
mtext("tau in (0, 1), xi = 1", side="3", outer=T)
dev.print(file="./plots/qscores/qscore-xipos-all.pdf", device=pdf)

tau <- seq(0.9, 1, 0.001)
xi <- 0.25
sig <- 1
quant <- sig / xi * ((1-tau)^(-xi) - 1)

par(mfrow=c(4, 4))
par(oma=c(2, 2, 2, 0))
par(mar=c(3, 3, 1, 1))
for (i in 1:n) {
	quant.score <- (quant - x[i]) * ((x[i] < quant) - tau)
	
	plot(tau, quant.score, type="l", main=paste("x = ", x[i]), xlab="", ylab="")

}
mtext("tau", side="1", outer=T)
mtext("quantile scores", side="2", outer=T)
mtext("tau in (0.90, 1), xi = 1", side="3", outer=T)
dev.print(file="./plots/qscores/qscore-xipos-highquant.pdf", device=pdf)

set.seed(2087)
xi <- 1
sig <- 1
n <- 16
x <- round(rgpd(n, 0, scale=sig, shape=xi), 3)
x <- sort(x)

tau <- seq(0, 1, 0.001)

quant <- sig / xi * ((1-tau)^(-xi) - 1)

par(mfrow=c(4, 4))
par(oma=c(2, 2, 2, 0))
par(mar=c(3, 3, 2, 1))
for (i in 1:n) {
	quant.score <- (quant - x[i]) * ((x[i] < quant) - tau)
	
	plot(tau, quant.score, type="l", main=paste("x = ", x[i]), xlab="", ylab="")

}
mtext("tau", side="1", outer=T)
mtext("quantile scores", side="2", outer=T)
mtext("tau in (0, 1), xi = 1", side="3", outer=T)
dev.print(file="./plots/qscores/qscore-xi1-all.pdf", device=pdf)

tau <- seq(0.9, 1, 0.001)
xi <- 1
sig <- 1
quant <- sig / xi * ((1-tau)^(-xi) - 1)

par(mfrow=c(4, 4))
par(oma=c(2, 2, 2, 0))
par(mar=c(3, 3, 1, 1))
for (i in 1:n) {
	quant.score <- (quant - x[i]) * ((x[i] < quant) - tau)
	
	plot(tau, quant.score, type="l", main=paste("x = ", x[i]), xlab="", ylab="")

}
mtext("tau", side="1", outer=T)
mtext("quantile scores", side="2", outer=T)
mtext("tau in (0.90, 1), xi = 1", side="3", outer=T)
dev.print(file="./plots/qscores/qscore-xi1-highquant.pdf", device=pdf)

set.seed(2087)
xi <- 1.5
sig <- 1
n <- 16
x <- round(rgpd(n, 0, scale=sig, shape=xi), 3)
x <- sort(x)

tau <- seq(0, 1, 0.001)

quant <- sig / xi * ((1-tau)^(-xi) - 1)

par(mfrow=c(4, 4))
par(oma=c(2, 2, 2, 0))
par(mar=c(3, 3, 2, 1))
for (i in 1:n) {
	quant.score <- (quant - x[i]) * ((x[i] < quant) - tau)
	
	plot(tau, quant.score, type="l", main=paste("x = ", x[i]), xlab="", ylab="")

}
mtext("tau", side="1", outer=T)
mtext("quantile scores", side="2", outer=T)
mtext("tau in (0, 1), xi = 1.5", side="3", outer=T)
dev.print(file="./plots/qscores/qscore-xi1.5-all.pdf", device=pdf)

tau <- seq(0.9, 1, 0.001)
xi <- 1.5
sig <- 1
quant <- sig / xi * ((1-tau)^(-xi) - 1)

par(mfrow=c(4, 4))
par(oma=c(2, 2, 2, 0))
par(mar=c(3, 3, 1, 1))
for (i in 1:n) {
	quant.score <- (quant - x[i]) * ((x[i] < quant) - tau)
	
	plot(tau, quant.score, type="l", main=paste("x = ", x[i]), xlab="", ylab="")

}
mtext("tau", side="1", outer=T)
mtext("quantile scores", side="2", outer=T)
mtext("tau in (0.90, 1), xi = 1.5", side="3", outer=T)
dev.print(file="./plots/qscores/qscore-xi1.5-highquant.pdf", device=pdf)

set.seed(2087)
xi <- -0.25
sig <- 1
n <- 16
x <- round(rgpd(n, 0, scale=sig, shape=xi), 3)
x <- sort(x)

tau <- seq(0, 1, 0.001)

quant <- sig / xi * ((1-tau)^(-xi) - 1)

par(mfrow=c(4, 4))
par(oma=c(2, 2, 2, 0))
par(mar=c(3, 3, 1, 1))
for (i in 1:n) {
	quant.score <- (quant - x[i]) * ((x[i] < quant) - tau)
	
	plot(tau, quant.score, type="l", main=paste("x = ", x[i]), xlab="", ylab="")

}
mtext("tau", side="1", outer=T)
mtext("quantile scores", side="2", outer=T)
mtext("tau in (0, 1), xi = -0.25", side="3", outer=T)
dev.print(file="./plots/qscores/qscore-xineg-all.pdf", device=pdf)

tau <- seq(0.9, 1, 0.001)
xi <- -0.25
sig <- 1
quant <- sig / xi * ((1-tau)^(-xi) - 1)

par(mfrow=c(4, 4))
par(oma=c(2, 2, 2, 0))
par(mar=c(3, 3, 1, 1))
for (i in 1:n) {
	quant.score <- (quant - x[i]) * ((x[i] < quant) - tau)
	
	plot(tau, quant.score, type="l", main=paste("x = ", x[i]), xlab="", ylab="")

}
mtext("tau", side="1", outer=T)
mtext("quantile scores", side="2", outer=T)
mtext("tau in (0.90, 1), xi = -0.25", side="3", outer=T)
dev.print(file="./plots/qscores/qscore-xineg-highquant.pdf", device=pdf)

set.seed(2087)
xi <- 1
sig <- 6
n <- 16
x <- round(rgpd(n, 0, scale=sig, shape=xi), 3)
x <- sort(x)

tau <- seq(0, 1, 0.001)

quant <- sig / xi * ((1-tau)^(-xi) - 1)

par(mfrow=c(4, 4))
par(oma=c(2, 2, 2, 0))
par(mar=c(3, 3, 2, 1))
for (i in 1:n) {
	quant.score <- (quant - x[i]) * ((x[i] < quant) - tau)
	
	plot(tau, quant.score, type="l", main=paste("x = ", x[i]), xlab="", ylab="")

}
mtext("tau", side="1", outer=T)
mtext("quantile scores", side="2", outer=T)
mtext("tau in (0, 1), xi = 1, sig=6", side="3", outer=T)
dev.print(file="./plots/qscores/qscore-x1sig6-all.pdf", device=pdf)

#### beta posteriors

load("thresh-1.RData")
par(mfrow = c(2, 2))
boxplot(beta.prob[[1]][1, , 1], main="betas in ex prob", outline=F)
abline(0, 0)
boxplot(beta.scale[[1]][1, , 1], main="betas in GDP scale", outline=F)
abline(0, 0)
boxplot(beta.shape[[1]][, 1, 1], main="betas in GPD shape", outline=F)  
abline(0, 0)