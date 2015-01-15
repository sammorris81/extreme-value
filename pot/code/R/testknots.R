nt <- 10
nknots <- 2
s.temp <- cbind(seq(0, 10, length=10), seq(0, 10, length=10))

set.seed(200)

phi <- 0
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10))
lines(knots.temp[2, , ], type="b", col="red")

par(mfrow=c(2, 3))
phi <- 0.50
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")

phi <- 0.70
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")

phi <- 0.80
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")

phi <- 0.90
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")

phi <- 0.95
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")

phi <- 0.99
knots.temp <- makeknotsTS(nt, nknots, s.temp, phi)
plot(knots.temp[1, , ], type="b", ylim=c(0, 10), xlim=c(0, 10),
     main=bquote(phi[w] == .(phi)), xlab="x", ylab="y")
lines(knots.temp[2, , ], type="b", col="red")