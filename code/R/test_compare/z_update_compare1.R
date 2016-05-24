library(inline)
library(SpatialTools)
library(fields)
library(microbenchmark)
#generate new z matrix
z1.R <- function(taug, tau, y, x.beta, mu, g, prec, lambda, zg) {
  nknots <- nrow(tau)
  nt     <- ncol(tau)
  z   <- matrix(NA, nknots, nt)
  mmm <- matrix(NA, nknots, nt)
  sd  <- matrix(NA, nknots, nt)
  for (t in 1:nt) {
    taug.t   <- sqrt(taug[, t])
    y.t      <- y[, t]
    x.beta.t <- x.beta[, t]
    mu.t     <- mu[, t]
    for (k in 1:nknots) {
      these <- which(g[, t] == k)

      r.1 <- (y.t[these] - x.beta.t[these])
      r.2 <- (y.t[-these] - mu.t[-these]) * taug.t[-these] / sqrt(tau[k, t])

      prec.11 <- prec[these, these, drop=F]
      prec.21 <- prec[-these, these, drop=F]

      mmmkt <- lambda * sum(r.1 %*% prec.11 + r.2 %*% prec.21)
      vvvkt <- 1 + lambda^2 * sum(prec.11)

      mmmkt <- mmmkt / vvvkt
      sdkt  <- 1 / sqrt(vvvkt * tau[k, t])
      sd[k, t] <- sdkt
      mmm[k, t] <- mmmkt
      z[k, t] <- abs(rnorm(1, mmm[k, t], sdkt))
      zg[these, t] <- z[k, t]
      mu.t <- x.beta.t + lambda * zg[, t]

    }
  }
  results <- list(z=z, mmm=mmm, sd=sd)
  return(results)
}

code1 <- '
  arma::mat tau_g_c = Rcpp::as<arma::mat>(taug);
  arma::mat tau_c = Rcpp::as<arma::mat>(tau);
  arma::mat y_c = Rcpp::as<arma::mat>(y);
  arma::mat xbeta_c = Rcpp::as<arma::mat>(x_beta);
  arma::mat mu_c = Rcpp::as<arma::mat>(mu);
  arma::mat g_c = Rcpp::as<arma::mat>(g);
  arma::mat prec_c = Rcpp::as<arma::mat>(prec);
  arma::mat prec_11; arma::mat prec_21;
  double lambda_c = Rcpp::as<double>(lambda);
  arma::mat zg_c = Rcpp::as<arma::mat>(zg);
  int nknots = tau_c.n_rows; int nt = tau_c.n_cols; int ns = tau_g_c.n_rows;
  int nparts;
  arma::mat z(nknots, nt);
  arma::uvec these; arma::uvec nthese; arma::vec tau_g_t (ns);
  arma::vec r_1; arma::vec r_2;
  arma::vec y_t(ns); arma::vec xbeta_t(ns); arma::vec mu_t(ns);
  arma::mat mmm(nknots, nt); arma::mat sd(nknots, nt);
  double mmmkt; double vvvkt; double sdkt; double zkt;

  for (int t = 0; t < nt; t++) {
    // extract the current values for day t
    tau_g_t = sqrt(tau_g_c.col(t));
    y_t = y_c.col(t); xbeta_t = xbeta_c.col(t); mu_t = mu_c.col(t);

    for (int k = 0; k < nknots; k++) {
      these  = find(g_c.col(t) == (k+1));
      nthese = find(g_c.col(t) != (k+1));  // nthese is "not these"
      nparts = these.n_elem;

      // recompute the mean for the day to account for updates to other knots
      mu_t = xbeta_t + lambda_c * zg_c.col(t);

      r_1 = (y_t(these) - xbeta_t(these));
      r_2 = (y_t(nthese) - mu_t(nthese)) % tau_g_t(nthese) / sqrt(tau_c(k, t));

      prec_11 = prec_c(these, these); prec_21 = prec_c(nthese, these);

      mmmkt = lambda_c * sum(r_1.t() * prec_11 + r_2.t() * prec_21);
      vvvkt = 1 + pow(lambda_c, 2) * accu(prec_11);

      mmmkt = mmmkt / vvvkt;
      sdkt  = 1 / sqrt(vvvkt * tau_c(k, t));
      sd(k, t) = sdkt;
      mmm(k, t) = mmmkt;
      zkt = (rnorm(1, mmmkt, sdkt))(0);  // vector of length 1
      zkt = std::abs(zkt);               // needs to have std::
      z(k, t) = zkt;

      // update the zg matrix
      for (int p = 0; p < nparts; p++) {
        zg_c(these(p), t) = zkt;
      }
    }

    mu_c.col(t) = mu_t;  // update mu
  }

  return Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("zg") = zg_c,
    Rcpp::Named("mmm") = mmm,
    Rcpp::Named("sd") = sd
  );
'
z.Rcpp1 <- cxxfunction(signature(taug="numeric", tau="numeric",
                                y="numeric", x_beta="numeric",
                                mu="numeric", g="numeric", prec="numeric",
                                lambda="numeric", zg="numeric"),
                                code1, plugin="RcppArmadillo")

set.seed(5)
ns <- 100
nknots <- 5
nt <- 10
s.temp <- cbind(runif(ns, 0, 10), runif(ns, 0, 10))
d.temp <- rdist(s.temp)
prec.temp <- simple.cov.sp(D=d.temp, sp.type="matern", sp.par=c(1, 1),
                           error.var=0, smoothness=0.5, finescale.var=0)
y <- matrix(rnorm(ns*nt, 0, 1), ns, nt)

x.beta <- matrix(rnorm(ns*nt, 0, 1), ns, nt)
g <- matrix(0, ns, nt)
for (t in 1:nt) {
  g[, t] <- sample(1:nknots, ns, replace=T)
}
lambda <- 2
lambda.1 <- 1
lambda.2 <- 1 / lambda^2
z <- matrix(abs(rnorm(nknots * nt)), nknots, nt)
zg <- matrix(0, ns, nt)
for (t in 1:nt) {
  zg[, t] <- z[g[, t], t]
}
tau <- matrix(exp(rnorm(nknots * nt)), nknots, nt)
taug <- matrix(0, ns, nt)
for (t in 1:nt) {
  taug[, t] <- tau[g[, t], t]
}
mu <- x.beta + lambda * zg

set.seed(10)
test.r <- z1.R(taug=taug, tau=tau, y=y, x.beta=x.beta, mu=mu, g=g,
              prec=prec.temp, lambda=lambda, zg=zg)

set.seed(10)
test.cpp <- z.Rcpp1(taug=taug, tau=tau, y=y, x_beta=x.beta, mu=mu, g=g,
                   prec=prec.temp, lambda=lambda, zg=zg)

print(test.cpp$sd / test.r$sd)
print(test.cpp$mmm / test.r$mmm)

microbenchmark(z1.R(taug=taug, tau=tau, y=y, x.beta=x.beta, mu=mu, g=g,
                   prec=prec.temp, lambda=lambda, zg=zg),
               z.Rcpp1(taug=taug, tau=tau, y=y, x_beta=x.beta, mu=mu, g=g,
                   prec=prec.temp, lambda=lambda, zg=zg),
               times=100)