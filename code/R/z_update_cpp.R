library(inline)
code <- '
  arma::mat tau_g_c = Rcpp::as<arma::mat>(taug);
  arma::mat tau_c = Rcpp::as<arma::mat>(tau);
  arma::mat y_c = Rcpp::as<arma::mat>(y);
  arma::mat xbeta_c = Rcpp::as<arma::mat>(x_beta);
  arma::mat mu_c = Rcpp::as<arma::mat>(mu);
  arma::mat g_c = Rcpp::as<arma::mat>(g);
  arma::mat prec_c = Rcpp::as<arma::mat>(prec);
  arma::mat prec_11; arma::mat prec_21;
  double lambda_1_c = Rcpp::as<double>(lambda_1);
  double lambda_2_c = Rcpp::as<double>(lambda_2);
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
      mu_t = xbeta_t + lambda_1_c * zg_c.col(t);

      r_1 = (y_t(these) - xbeta_t(these));
      r_2 = (y_t(nthese) - mu_t(nthese)) % tau_g_t(nthese) / sqrt(tau_c(k, t));

      prec_11 = prec_c(these, these); prec_21 = prec_c(nthese, these);

      mmmkt = lambda_1_c * sum(r_1.t() * prec_11 + r_2.t() * prec_21);
      vvvkt = lambda_2_c + accu(prec_11);

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
z.Rcpp <- cxxfunction(signature(taug="numeric", tau="numeric",
                                y="numeric", x_beta="numeric",
                                mu="numeric", g="numeric", prec="numeric",
                                lambda_1="numeric", lambda_2="numeric",
                                zg="numeric"),
                                code, plugin="RcppArmadillo")

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
      zkt = std::abs(zkt);                      // needs to have std::
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