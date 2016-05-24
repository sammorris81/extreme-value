library(fields)
library(colorspace)
s1 <- s2 <- seq(0.001, 0.999, length = 150)
s <- as.matrix(expand.grid(s1, s2))

nknots <- 8
w1 <- runif(8, 0, 1)
w2 <- runif(8, 0, 1)

w <- cbind(w1, w2)

s.star <- qnorm(s)
w.star <- qnorm(w)

d.star <- rdist(s.star, w.star)
d <- rdist(s, w)

get.min <- function(vec) {
  return (which(vec == min(vec)))
}

part.star <- apply(d.star, 1, get.min)
part      <- apply(d, 1, get.min)

plot(s)
plot(s, pch = 19, col = rainbow_hcl(8)[part])
points(w, pch = 20, col = 1)

par(mfrow = c(1, 2))
plot(s.star, pch = 19, col = rainbow_hcl(8)[part])
plot(s.star, pch = 19, col = rainbow_hcl(8)[part.star])

