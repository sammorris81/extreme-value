nt <- 10
nknots <- 2
s.temp <- cbind(seq(0, 10, length=10), seq(0, 10, length=10))

set.seed(10)
phi <- 0
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10))
lines(knots.temp[2, , ], type="b", col="red")

par(mfrow=c(2, 3))
set.seed(20)
phi <- 0.50
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")

set.seed(30)
phi <- 0.70
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")

set.seed(40)
phi <- 0.80
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")

set.seed(50)
phi <- 0.90
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")

set.seed(60)
phi <- 0.95
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")

set.seed(70)
phi <- 0.99
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")

# making TS for tau
nt <- 1000
nknots <- 2
tau.alpha <- 3
tau.beta <- 3
xplot <- seq(0, 5, 0.01)
yplot <- dgamma(xplot, tau.alpha, tau.beta)
par(mfrow=c(2, 2))

set.seed(10)
phi <- 0
tau.temp <- maketauTS(nt, nknots, tau.alpha, tau.beta, phi)
for (i in 1:2) {
  y.max <- max(hist(tau.temp[i, ], plot=F)$density, yplot)
  hist(tau.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[tau] == .(phi)),
       xlab=bquote(tau))
  lines(xplot, yplot)
  plot(tau.temp[i, ], type="l", ylab=bquote(tau))
}

set.seed(20)
phi <- 0.50
tau.temp <- maketauTS(nt, nknots, tau.alpha, tau.beta, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(tau.temp[i, ], plot=F)$density, yplot)
  hist(tau.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[tau] == .(phi)),
       xlab=bquote(tau))
  lines(xplot, yplot)
  plot(tau.temp[i, ], type="l", ylab=bquote(tau))
}

set.seed(30)
phi <- 0.70
tau.temp <- maketauTS(nt, nknots, tau.alpha, tau.beta, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(tau.temp[i, ], plot=F)$density, yplot)
  hist(tau.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[tau] == .(phi)),
       xlab=bquote(tau))
  lines(xplot, yplot)
  plot(tau.temp[i, ], type="l", ylab=bquote(tau))
}

set.seed(40)
phi <- 0.80
tau.temp <- maketauTS(nt, nknots, tau.alpha, tau.beta, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(tau.temp[i, ], plot=F)$density, yplot)
  hist(tau.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[tau] == .(phi)),
       xlab=bquote(tau))
  lines(xplot, yplot)
  plot(tau.temp[i, ], type="l", ylab=bquote(tau))
}

set.seed(50)
phi <- 0.90
tau.temp <- maketauTS(nt, nknots, tau.alpha, tau.beta, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(tau.temp[i, ], plot=F)$density, yplot)
  hist(tau.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[tau] == .(phi)),
       xlab=bquote(tau))
  lines(xplot, yplot)
  plot(tau.temp[i, ], type="l", ylab=bquote(tau))
}

set.seed(60)
phi <- 0.95
tau.temp <- maketauTS(nt, nknots, tau.alpha, tau.beta, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(tau.temp[i, ], plot=F)$density, yplot)
  hist(tau.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[tau] == .(phi)),
       xlab=bquote(tau))
  lines(xplot, yplot)
  plot(tau.temp[i, ], type="l", ylab=bquote(tau))
}

set.seed(70)
phi <- 0.99
tau.temp <- maketauTS(nt, nknots, tau.alpha, tau.beta, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(tau.temp[i, ], plot=F)$density, yplot)
  hist(tau.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[tau] == .(phi)),
       xlab=bquote(tau))
  lines(xplot, yplot)
  plot(tau.temp[i, ], type="l", ylab=bquote(tau))
}

# making TS for z
nt <- 500
nknots <- 2
xplot <- seq(0, 10, 0.01)
tau.scalar <- 2
yplot <- dnorm(xplot, 0, sqrt(1 / tau.scalar)) / 0.5
tau <- matrix(tau, nrow=nknots, ncol=nt)  # for the function to make the z

set.seed(10)
phi <- 0
z.temp <- makezTS(nt, nknots, tau, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(z.temp[i, ], plot=F)$density, yplot)
  hist(z.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[z] == .(phi)),
       xlab="z")
  lines(xplot, yplot)
  plot(z.temp[i, ], type="l", ylab="z")
}

set.seed(20)
phi <- 0.50
z.temp <- makezTS(nt, nknots, tau, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(z.temp[i, ], plot=F)$density, yplot)
  hist(z.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[z] == .(phi)),
       xlab="z")
  lines(xplot, yplot)
  plot(z.temp[i, ], type="l", ylab="z")
}

set.seed(30)
phi <- 0.70
z.temp <- makezTS(nt, nknots, tau, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(z.temp[i, ], plot=F)$density, yplot)
  hist(z.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[z] == .(phi)),
       xlab="z")
  lines(xplot, yplot)
  plot(z.temp[i, ], type="l", ylab="z")
}

set.seed(40)
phi <- 0.80
z.temp <- makezTS(nt, nknots, tau, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(z.temp[i, ], plot=F)$density, yplot)
  hist(z.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[z] == .(phi)),
       xlab="z")
  lines(xplot, yplot)
  plot(z.temp[i, ], type="l", ylab="z")
}

set.seed(50)
phi <- 0.90
z.temp <- makezTS(nt, nknots, tau, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(z.temp[i, ], plot=F)$density, yplot)
  hist(z.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[z] == .(phi)),
       xlab="z")
  lines(xplot, yplot)
  plot(z.temp[i, ], type="l", ylab="z")
}

set.seed(60)
phi <- 0.95
z.temp <- makezTS(nt, nknots, tau, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(z.temp[i, ], plot=F)$density, yplot)
  hist(z.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[z] == .(phi)),
       xlab="z")
  lines(xplot, yplot)
  plot(z.temp[i, ], type="l", ylab="z")
}

set.seed(70)
phi <- 0.99
z.temp <- makezTS(nt, nknots, tau, phi)
par(mfrow=c(2, 2))
for (i in 1:2) {
  y.max <- max(hist(z.temp[i, ], plot=F)$density, yplot)
  hist(z.temp[i, ], freq=F, ylim=c(0, y.max), main=bquote(phi[z] == .(phi)),
       xlab="z")
  lines(xplot, yplot)
  plot(z.temp[i, ], type="l", ylab="z")
}

