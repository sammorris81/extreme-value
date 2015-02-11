updateLambda1 <- function(x.beta, zg, y, prec, taug) {
  lll <- c(0, 0)
  mmm <- c(-1, 1)
  res.l <- y - (x.beta + mmm * zg)
  for (l in 1:length(lll)) {
    lll[l] <- -0.5 * sum(rss(prec, y=sqrt(taug) * res.l[l]))
  }
  lambda.1 <- sample(mmm, 1, prob=exp(lll - max(lll)))
}

updateLambda2 <- function(lambda.a, lambda.b, z, tau) {
  nknots <- nrow(z)
  nt     <- ncol(z)
  aaa    <- lambda.a + 0.5 * nknots * nt
  bbb    <- lambda.b + 0.5 * sum(z^2 * tau)
  lambda.2 <- rgamma(1, aaa, bbb)
}