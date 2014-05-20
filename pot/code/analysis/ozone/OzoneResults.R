rm(list=ls())

load("OzoneMCMCnothresh.RData")
load("OzoneMCMCthresh.RData")

delta <- matrix(0, nrow=5000, ncol=4)
delta[, 1] <- fit.normal.nothresh$delta[15001:20000]
delta[, 2] <- fit.normal.thresh90$delta[25001:30000]
delta[, 3] <- fit.skew.nothresh$delta[15001:20000]
delta[, 4] <- fit.skew.thresh90$delta[25001:30000]
colnames(delta) <- c("Normal (T=0)", "Normal (T=0.9)", "Skew (T=0)", "Skew (T=0.9)")
boxplot(delta)

alpha <- matrix(0, nrow=5000, ncol=4)
alpha[, 1] <- fit.normal.nothresh$alpha[15001:20000]
alpha[, 2] <- fit.normal.thresh90$alpha[25001:30000]
alpha[, 3] <- fit.skew.nothresh$alpha[15001:20000]
alpha[, 4] <- fit.skew.thresh90$alpha[25001:30000]
colnames(alpha) <- c("Normal (T=0)", "Normal (T=0.9)", "Skew (T=0)", "Skew (T=0.9)")
boxplot(alpha)


