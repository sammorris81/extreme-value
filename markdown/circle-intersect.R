library(plotrix)
plot(0, 0, ylim = c(-1.5, 1.5), xlim = c(-1.5, 1.5))

xplot <- seq(-1, 1, length = 1000)
yplot <- sqrt(1 - xplot^2)
lines(xplot, yplot)

yplot <- -sqrt(1 - xplot^2)
lines(xplot, yplot)

v1 <- c(-.2, .7)
v2 <- c(-.3, .9)

points(1, 0, pch = 16)
points(-1, 0, pch = 16)
points(v1[1], v1[2])
points(v2[1], v2[2])

midpoint.x <- (v1[1] + v2[1]) / 2
midpoint.y <- (v1[2] + v2[2]) / 2

slope <- -(v1[1] - v2[1]) / (v1[2] - v2[2])
int <- slope * -midpoint.x + midpoint.y
xplot <- seq(-1.5, 1.5)
yplot <- int + slope * xplot
lines(xplot, yplot)

# circle: x^2 + y^2 = 1

# y coordinates
y1 <- (2 * int + sqrt(4 * int^2 - 4 * (1 + slope^2) * (int^2 - slope^2))) / (2 + 2 * slope^2)
y2 <- (2 * int - sqrt(4 * int^2 - 4 * (1 + slope^2) * (int^2 - slope^2))) / (2 + 2 * slope^2)

x1 <- (y1 - int) / slope
x2 <- (y2 - int) / slope

points(x1, y1)
points(x2, y2)



# look at overlapping circles
xplot <- seq(-1, 1, length = 10000)
yplot.pos <- sqrt(1 - xplot^2)
yplot.neg <- -sqrt(1 - xplot^2)

s1 <- c(-1, 0)
s2 <- c(1, 0)
points(s1[1], s1[2], pch = 16)
points(s2[1], s2[2], pch = 16)

v <- c(-0.4, 0.6)
d1 <- sqrt((s1[1] - v[1])^2 + (s1[2] - v[2])^2)

d2 <- sqrt((s2[1] - v[1])^2 + (s2[2] - v[2])^2)

xplot1 <- xplot * d1 + s1[1]
yplot1.pos <- yplot.pos * d1 + s1[2]
yplot1.neg <- yplot.neg * d1 + s1[2]

xplot2 <- xplot * d2 + s2[1]
yplot2.pos <- yplot.pos * d2 + s2[2]
yplot2.neg <- yplot.neg * d2 + s2[2]

ylims <- range(c(yplot.pos, yplot1.pos, yplot2.pos, yplot.neg, yplot1.neg, yplot2.neg))
xlims <- range(c(xplot, xplot1, xplot2))

plot(0, 0, ylim = ylims, xlim = xlims, type = "n")
abline(h = 0)
lines(xplot, yplot.pos)
lines(xplot, yplot.neg)
lines(xplot1, yplot1.pos)
lines(xplot1, yplot1.neg)
lines(xplot2, yplot2.pos)
lines(xplot2, yplot2.neg)


# look at overlapping circles
xplot <- seq(-1, 1, length = 10000)
yplot.pos <- sqrt(1 - xplot^2)
yplot.neg <- -sqrt(1 - xplot^2)

s1 <- c(-1, 0)
s2 <- c(1, 0)
points(s1[1], s1[2], pch = 16)
points(s2[1], s2[2], pch = 16)

v <- c(0, 0)
d1 <- sqrt((s1[1] - v[1])^2 + (s1[2] - v[2])^2)

d2 <- sqrt((s2[1] - v[1])^2 + (s2[2] - v[2])^2)

xplot1 <- xplot * d1 + s1[1]
yplot1.pos <- yplot.pos * d1 + s1[2]
yplot1.neg <- yplot.neg * d1 + s1[2]

xplot2 <- xplot * d2 + s2[1]
yplot2.pos <- yplot.pos * d2 + s2[2]
yplot2.neg <- yplot.neg * d2 + s2[2]

ylims <- range(c(yplot.pos, yplot1.pos, yplot2.pos, yplot.neg, yplot1.neg, yplot2.neg))
xlims <- range(c(xplot, xplot1, xplot2))

plot(0, 0, ylim = ylims, xlim = xlims, type = "n")
abline(h = 0)
lines(xplot, yplot.pos)
lines(xplot, yplot.neg)
lines(xplot1, yplot1.pos)
lines(xplot1, yplot1.neg)
lines(xplot2, yplot2.pos)
lines(xplot2, yplot2.neg)

reps <- 100000
area <- 0
for (i in 1:reps) {
  r <- runif(1, 0, 1)
  theta <- runif(1, 0, 2 * pi)
  x <- sqrt(r) * cos(theta)
  y <- sqrt(r) * sin(theta)
  if (x < 0) {
    wedge <- ((x + 1)^2 + y^2) >= 1
    area <- area + wedge / reps
    if (wedge) {
      points(x, y, pch = 19)
    }
  } else if (x >= 0) {
    wedge <- ((x - 1)^2 + y^2) >= 1
    area <- area + wedge / reps
    if (wedge) {
      points(x, y, pch = 19)
    }
  }
}
area * pi