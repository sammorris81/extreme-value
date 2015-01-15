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
  double lambda_c = Rcpp::as<double>(lambda);
  arma::mat zg_c = Rcpp::as<arma::mat>(zg);
  int nknots = tau_c.n_rows; int nt = tau_c.n_cols; int ns = tau_g_c.n_rows;
  int nparts;
  arma::mat z(nknots, nt);
  arma::uvec these; arma::uvec nthese; arma::vec tau_g_t (ns);
  arma::vec r_1; arma::vec r_2;
  arma::vec y_t(ns); arma::vec xbeta_t(ns); arma::vec mu_t(ns);
  arma::mat mmm(nknots, nt); arma::mat vvv(nknots, nt);
  double mmmkt; double vvvkt; double zkt;
  
  for (int t = 0; t < nt; t++) {
    tau_g_t = sqrt(tau_g_c.col(t));
    y_t = y_c.col(t); xbeta_t = xbeta_c.col(t); mu_t = mu_c.col(t);
    for (int k = 0; k < nknots; k++) {
      these = find(g_c.col(t) == (k+1));
      nparts = these.n_elem;
      nthese = find(g_c.col(t) != (k+1));
      
      r_1 = (y_t(these) - xbeta_t(these)) % tau_g_t(these);
      r_2 = (y_t(nthese) - mu_t(nthese)) % tau_g_t(nthese);
      
      prec_11 = prec_c(these, these); prec_21 = prec_c(nthese, these);

      mmmkt = lambda_c * sqrt(tau_c(k, t)) * sum(r_1.t() * prec_11 + r_2.t() * prec_21);
      vvvkt = tau_c(k, t) + pow(lambda_c, 2) * tau_c(k, t) * accu(prec_11);
      
      vvvkt = 1 / vvvkt;
      vvv(k, t) = vvvkt;
      mmmkt = vvvkt * mmmkt;
      mmm(k, t) = mmmkt;
      zkt = (rnorm(1, mmmkt, sqrt(vvvkt)))(0);  // vector of length 1
      zkt = std::abs(zkt);                      // std is needed to give absolute value
      z(k, t) = zkt;
      
      for (int p = 0; p < nparts; p++) {
        zg_c(these(p), t) = zkt;
      }
      mu_t = xbeta_t + lambda_c * zg_c.col(t);
      
    }
    
    mu_c.col(t) = mu_t;
  }

  return Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("zg") = zg_c,
    Rcpp::Named("mmm") = mmm,
    Rcpp::Named("vvv") = vvv
  );
'
z.Rcpp <- cxxfunction(signature(taug="numeric", tau="numeric",
                                y="numeric", x_beta="numeric",
                                mu="numeric", g="numeric", prec="numeric",
                                lambda="numeric", zg="numeric"), 
                      code, plugin="RcppArmadillo")