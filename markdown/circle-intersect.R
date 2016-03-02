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
