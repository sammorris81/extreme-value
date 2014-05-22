rm(list=ls())

load("OzoneMCMCnothresh.RData")
load("OzoneMCMCthresh.RData")

start <- 25001
end   <- 30000

# posterior distribution delta
delta <- matrix(0, nrow=5000, ncol=4)
delta[, 1] <- fit.normal.nothresh$delta[start:end]
delta[, 2] <- fit.normal.thresh90$delta[start:end]
delta[, 3] <- fit.skew.nothresh$delta[start:end]
delta[, 4] <- fit.skew.thresh90$delta[start:end]
colnames(delta) <- c("Normal (T=0)", "Normal (T=0.9)", "Skew (T=0)", "Skew (T=0.9)")
boxplot(delta, main="Delta posterior")

# iteration plot delta
par(mfrow=c(2, 1))
plot(fit.skew.nothresh$delta, type="l", main="Delta iteration plot (T=0)", ylab="delta")
abline(v=25000, lty=2)
plot(fit.skew.thresh90$delta, type="l", main="Delta iteration plot (T=0.9)", ylab="delta")
abline(v=25000, lty=2)

# posterior distribution alpha
alpha <- matrix(0, nrow=5000, ncol=4)
alpha[, 1] <- fit.normal.nothresh$alpha[start:end]
alpha[, 2] <- fit.normal.thresh90$alpha[start:end]
alpha[, 3] <- fit.skew.nothresh$alpha[start:end]
alpha[, 4] <- fit.skew.thresh90$alpha[start:end]
colnames(alpha) <- c("Normal (T=0)", "Normal (T=0.9)", "Skew (T=0)", "Skew (T=0.9)")
boxplot(alpha, main="Alpha posterior")

# iteration plot alpha
par(mfrow=c(2, 1))
plot(fit.skew.nothresh$alpha, type="l", main="Alpha iteration plot (Skew, T=0)", ylab="alpha")
abline(v=25000, lty=2)
plot(fit.skew.thresh90$alpha, type="l", main="Alpha iteration plot (Skew, T=0.9)", ylab="alpha")
abline(v=25000, lty=2)
plot(fit.normal.nothresh$alpha, type="l", main="Alpha iteration plot (Normal, T=0)", ylab="alpha")
abline(v=25000, lty=2)
plot(fit.normal.thresh90$alpha, type="l", main="Alpha iteration plot (Normal, T=0.9)", ylab="alpha")
abline(v=25000, lty=2)

colnames(fit.normal.nothresh$beta) <- colnames(fit.normal.thresh90$beta) <- colnames(fit.skew.nothresh$beta) <- colnames(fit.skew.thresh90$beta) <- c("Int", "Long", "Lat", "CMAQ")

par(mfrow=c(2, 2))

boxplot(fit.normal.nothresh$beta[start:end, ], main="Beta posterior: Normal (T=0)")
boxplot(fit.normal.thresh90$beta[start:end, ], main="Beta posterior: Normal (T=0.9)")
boxplot(fit.skew.nothresh$beta[start:end, ], main="Beta posterior: Skew (T=0)")
boxplot(fit.skew.thresh90$beta[start:end, ], main="Beta posterior: Skew (T=0.9)")

beta.int <- matrix(0, nrow=5000, ncol=4)
beta.int[, 1] <- fit.normal.nothresh$beta[start:end, 1]
beta.int[, 2] <- fit.normal.thresh90$beta[start:end, 1]
beta.int[, 3] <- fit.skew.nothresh$beta[start:end, 1]
beta.int[, 4] <- fit.skew.thresh90$beta[start:end, 1]
colnames(beta.int) <- c("Normal (T=0)", "Normal (T=0.9)", "Skew (T=0)", "Skew (T=0.9)")
boxplot(beta.int, main="Beta intercept posterior")

beta.int <- matrix(0, nrow=5000, ncol=4)
beta.int[, 1] <- fit.normal.nothresh$beta[start:end, 2]
beta.int[, 2] <- fit.normal.thresh90$beta[start:end, 2]
beta.int[, 3] <- fit.skew.nothresh$beta[start:end, 2]
beta.int[, 4] <- fit.skew.thresh90$beta[start:end, 2]
colnames(beta.int) <- c("Normal (T=0)", "Normal (T=0.9)", "Skew (T=0)", "Skew (T=0.9)")
boxplot(beta.int, main="Beta Longitude posterior")

beta.int <- matrix(0, nrow=5000, ncol=4)
beta.int[, 1] <- fit.normal.nothresh$beta[start:end, 3]
beta.int[, 2] <- fit.normal.thresh90$beta[start:end, 3]
beta.int[, 3] <- fit.skew.nothresh$beta[start:end, 3]
beta.int[, 4] <- fit.skew.thresh90$beta[start:end, 3]
colnames(beta.int) <- c("Normal (T=0)", "Normal (T=0.9)", "Skew (T=0)", "Skew (T=0.9)")
boxplot(beta.int, main="Beta Latitude posterior")

beta.int <- matrix(0, nrow=5000, ncol=4)
beta.int[, 1] <- fit.normal.nothresh$beta[start:end, 4]
beta.int[, 2] <- fit.normal.thresh90$beta[start:end, 4]
beta.int[, 3] <- fit.skew.nothresh$beta[start:end, 4]
beta.int[, 4] <- fit.skew.thresh90$beta[start:end, 4]
colnames(beta.int) <- c("Normal (T=0)", "Normal (T=0.9)", "Skew (T=0)", "Skew (T=0.9)")
boxplot(beta.int, main="Beta CMAQ posterior")
