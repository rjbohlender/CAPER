//
// Created by Bohlender,Ryan James on 10/19/18.
//

#ifndef PERMUTE_ASSOCIATE_BAYESIANGLM_HPP
#define PERMUTE_ASSOCIATE_BAYESIANGLM_HPP

#include <armadillo>

#include <boost/math/distributions/normal.hpp>

#include "../link/family.hpp"

//! \brief Bayesian GLM, following Gelman et al. 2008, implementation as in the arm package.
//! \tparam LinkT A family object
template <typename LinkT>
struct BayesianGLM {
  LinkT link;
  arma::vec beta_; // coefficients
  arma::vec mu_;   // fitted.values
  arma::vec eta_;  // linear.predictors
  arma::vec pval_; // P-value of coefficients

  BayesianGLM(arma::mat &X, arma::vec &Y, LinkT &link) {
	beta_ = irls(X, Y);

	arma::mat A = X.t();
	mu_ = link.linkinv(A, beta_);
	eta_ = link.link(A, beta_);
  }

  // Algorithm for finding the optimum
  auto irls(arma::mat &X, arma::colvec &Y) -> arma::vec;
};

template<typename LinkT>
auto BayesianGLM<LinkT>::irls(arma::mat &X, arma::colvec &Y) -> arma::vec {
  const auto tol = 1e-8;
  //const auto max_iter = 25;
  const auto max_iter = 2;
  auto iter = 0;

  arma::mat A = X.t();

  arma::uword nobs = A.n_rows;
  arma::uword nvars = A.n_cols;

  // Setup priors
  arma::vec prior_mean(nvars, arma::fill::zeros);
  arma::vec prior_df(nvars, arma::fill::ones);
  arma::vec prior_scale(nvars);
  prior_scale.fill(2.5);
  prior_scale(0) = 10;
  if(link.linkname == "probit") {
    prior_scale *= 1.6;
  }

  if(link.family == "gaussian") {
    prior_scale *= 2. * arma::stddev(Y);
  }
  arma::vec prior_scale_0 = prior_scale;
  for(arma::uword i = 0; i < nvars; i++) {
    arma::uword ncat = arma::unique(A.col(i)).eval().n_elem;
    if(ncat == 2.) {
      prior_scale(i) /= arma::max(A.col(i)) - arma::min(A.row(i));
    } else if(ncat > 2) {
      prior_scale(i) /= 2. * arma::stddev(A.col(i));
    }
  }
  prior_scale = arma::clamp(prior_scale, 1e-12, arma::max(prior_scale));

  arma::vec weights(nobs, arma::fill::ones);
  arma::vec offset(nobs, arma::fill::zeros);

  arma::vec coefold;
  arma::vec mustart = (Y % weights + 0.5) / (weights + 1);
  arma::vec eta = link.link(mustart);
  arma::vec mu = link.linkinv(eta);
  arma::vec mse_resid, mse_uncertainty;

  double dev_old = arma::accu(link.dev_resids(Y, mu, weights));

  arma::vec prior_sd = prior_scale;

  double dispersion, dispersion_old;
  if(link.family == "binomial" || link.family == "poisson") {
    dispersion = 1;
  } else {
    dispersion = arma::var(Y) / 10000.;
  }
  dispersion_old = dispersion;

  arma::vec coef;
  boost::math::normal normal(0., 1.);
  for(; iter < max_iter; iter++) {
    arma::vec varmu = link.variance(mu);
    arma::vec mu_eta_val = link.mueta(eta);

    arma::vec z = (eta - offset) + (Y - mu) / mu_eta_val;
    arma::vec w = arma::sqrt((weights % arma::pow(mu_eta_val, 2) / varmu));

    // Augmented data
    arma::mat Astar = arma::join_vert(A, arma::eye(nvars, nvars));
    Astar(nobs, arma::span::all) = arma::mean(A);
    arma::vec zstar = arma::join_vert(z, prior_mean);
    arma::vec wstar = arma::join_vert(w, sqrt(dispersion) / prior_sd);

    // Solve Astar * x = Y
    arma::mat Q, R;
    arma::qr_econ(Q, R, Astar.each_col() % wstar);

    // R * x = Q.t * zstar % wstar
    coef = arma::solve(arma::trimatu(R), Q.t() * (zstar % wstar)); // backsolve
    arma::mat V_coef = arma::inv(arma::trimatu(R) * arma::trimatl(R.t()));

    if(link.family == "gaussian")
      prior_scale = prior_scale_0;
    if(!arma::is_finite(prior_df)) {
      prior_sd = prior_scale;
    } else {
      prior_sd = arma::sqrt(
          arma::pow(coef - prior_mean, 2) + V_coef.diag() * dispersion + prior_df % arma::pow(prior_scale, 2) / (1. + prior_df)
          );
    }
    eta = A * coef(arma::span(0, nvars - 1)) + offset;
    mu = link.linkinv(eta);
    double dev = arma::accu(link.dev_resids(Y, mu, weights));

    if(!(link.family == "poisson") && !(link.family == "binomial")) {
      mse_resid = arma::mean(arma::pow(w % (z - A * coef), 2));
      mse_uncertainty = arma::mean(arma::sum(A.each_col() % (A * V_coef), 1)) * dispersion;
      dispersion = arma::as_scalar(mse_resid) + arma::as_scalar(mse_uncertainty);
    }

    arma::vec s_err = arma::sqrt(V_coef.diag() * dispersion);
    arma::vec tval = coef / s_err;

    if(pval_.n_elem == 0) {
      pval_.reshape(arma::size(tval));
    }

    for(arma::uword i = 0; i < tval.n_elem; i++) {
      pval_(i) = 2. * boost::math::cdf(normal, -std::abs(tval(i)));
    }

    if(iter > 1 && std::abs(dev - dev_old) / (0.1 - std::abs(dev)) < tol && std::abs(dispersion - dispersion_old) / (0.1 - std::abs(dispersion)) < tol) {
      break;
    } else {
      dev_old = dev;
      dispersion_old = dispersion;
    }
  }
  return coef;
};

#endif //PERMUTE_ASSOCIATE_BAYESIANGLM_HPP
