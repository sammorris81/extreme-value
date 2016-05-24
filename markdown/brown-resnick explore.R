x <- seq(0, 10, length = 100)

set.seed(2000)

data <- matrix(NA, 5, length(x))

# as smooth -> 2, smooths out like smith
data[1, ] <- rmaxstab(1, x, cov.mod = "brown", range = 3, smooth = 1.99)

data[2, ] <- rmaxstab(1, x, cov.mod = "brown", range = 3, smooth = 1.5)

data[3, ] <- rmaxstab(1, x, cov.mod = "brown", range = 3, smooth = 1)

data[4, ] <- rmaxstab(1, x, cov.mod = "brown", range = 3, smooth = 0.5)

data[5, ] <- rmaxstab(1, x, cov.mod = "brown", range = 3, smooth = 0.25)

matplot(t(data[-4, ]), type = "b")
matplot(t(data[4:5, ]))
plot(x, data[4, ], type = "l")
# semivariogram gamma(h) = (|h| / range)^smooth