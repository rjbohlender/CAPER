//
// Created by Bohlender,Ryan James on 8/2/18.
//

#include <cassert>
#include <iostream>
#include "logistic_regression.hpp"
#include <boost/math/distributions/chi_squared.hpp>

/// Armadillo stores matrices in column major format. Iteration proceeds along columns.

LogisticRegression::LogisticRegression(arma::mat &Xmat, arma::colvec &Yvec)
: theta_(arma::rowvec(Xmat.n_rows, arma::fill::ones)) {
  // Switch to gradient descent for large samples
  try{
    theta_ = irls_svdnewton(Xmat, Yvec);
  } catch(std::exception &e) {
    std::cerr << "IRLS failed; Using gradient descent." << std::endl;
    theta_ = gradient_descent(Xmat, Yvec);
  }

  calculate_probability(Xmat);
  calculate_odds(Xmat);
  calculate_mean(Xmat);
  calculate_eta();
}

auto LogisticRegression::h(arma::mat &Xmat, arma::rowvec &t) -> arma::rowvec {
  return 1. / (1. + arma::exp(-(t * Xmat)));
}

//! \brief The log likelihood function.
//! \param Xmat Fixed effects, or combined fixed and random effects
//! \param Yvec Response vector
//! \param t Effect of each variable
//! \return
auto LogisticRegression::cost(arma::mat &Xmat, arma::colvec &Yvec, arma::rowvec &t) -> double {
  double m = Yvec.n_rows;

  arma::mat ret = -1. / m * (arma::log(h(Xmat, t)) * Yvec + arma::log(1 - h(Xmat, t)) * (1 - Yvec));

  assert(ret.n_rows == 1);
  assert(ret.n_cols == 1);

  return arma::as_scalar(ret);
}

auto LogisticRegression::gradient_descent(arma::mat &Xmat, arma::colvec &Yvec) -> arma::rowvec {
  std::cerr << "Running gradient descent.\n";
  auto iterations = 0ull;
  auto max_iter = 1000000ull;
  auto alpha = 0.005; // Learning rate
  auto tol = 1e-10;
  auto m = static_cast<double>(Xmat.n_cols);
  auto t = arma::rowvec(Xmat.n_rows, arma::fill::zeros);
  auto grad = t;
  do {
    // Vectorized update
    grad = alpha * (Xmat * (h(Xmat, theta_).t() - Yvec)).t() / m;
    t -= grad;


    iterations++;
  } while(iterations < max_iter && arma::any(grad > tol));
  return t;
}

auto LogisticRegression::calculate_odds(arma::mat &Xmat) -> void {
  if(arma::any(mu_ <= 0)) {
    std::cerr << "mu_ contains zeros\n";
  }
  odds_ = mu_ / (1. - mu_);
}

auto LogisticRegression::get_odds() -> arma::vec {
  return odds_;
}

auto LogisticRegression::calculate_mean(arma::mat &Xmat) -> void {
  mean_  = arma::mean(mu_);
}

auto LogisticRegression::get_mean() -> double {
  return mean_;
}

auto LogisticRegression::get_probability() -> arma::vec& {
  return mu_;
}

auto LogisticRegression::calculate_probability(arma::mat &Xmat) -> void {
  // Prevent numerical issues
  mu_ = arma::clamp(h(Xmat, theta_).t(), 0.00001, 0.99999);
}

auto LogisticRegression::get_eta() -> arma::vec& {
  return eta_;
}

auto LogisticRegression::calculate_eta() -> void {
  eta_ = arma::log(mu_ / (1. - mu_));
}

//! \brief Performs Iteratively Reweighted Least Squares
//! \param Xmat Features
//! \param Yvec Response vector
//!
//! Algorithm from: https://bwlewis.github.io/GLM/#svdnewton
//!
auto LogisticRegression::irls_svdnewton(arma::mat &Xmat, arma::colvec &Yvec) -> arma::rowvec {
  // Iteration params
  std::cerr << "Running IRLS.\n";
  const auto tol = 1e-8;
  const auto max_iter = 25;
  auto update = 0.;
  auto iter = 0;

  // Aliases
  arma::mat A = Xmat.t();
  arma::vec b = Yvec;

  arma::uword m = A.n_rows;
  arma::uword n = A.n_cols;

  // SVD
  arma::mat U, V;
  arma::vec S;

  arma::svd_econ(U, S, V, A);

  // Matrices and Vectors
  arma::vec t(m, arma::fill::zeros);
  arma::vec s(n, arma::fill::zeros);
  arma::vec s_old;
  arma::vec weights(m, arma::fill::ones);

  // Link Function -- Inverse logit and derivative
  auto g = [](arma::vec x) -> arma::vec {
    return 1. / (1. + arma::exp(-x));
  };
  auto gprime = [](arma::vec x) -> arma::vec {
    return arma::exp(x) / arma::pow((arma::exp(x) + 1.), 2);
  };

  do {
    arma::vec gv = g(t);
    arma::vec varg = gv % (1. - gv);
    arma::vec gvp = gprime(t);

    arma::vec z(m, arma::fill::zeros);
    arma::vec W(m, arma::fill::zeros);

    z = t + (b - gv) / gvp;
    W = weights % arma::pow(gvp, 2) / varg;


    s_old = s;

    arma::mat C;
    bool success = arma::chol(C, U.t() * (U.each_col() % W));
    if(!success) {
      std::cerr << z.t();
      std::cerr << W.t();
      std::cerr << gv.t();
      std::cerr << varg.t();
      std::cerr << gvp.t();
    }

    s = solve(arma::trimatl(C.t()), U.t() * (W % z));
    s = solve(arma::trimatu(C), s);

    t = U * s;

    update = arma::norm(s - s_old);

    iter++;
  } while(iter < max_iter && update > tol);

  return (V * (arma::diagmat(1. / S) * (U.t() * t))).t();
}

auto LogisticRegression::hessian(arma::mat &X) -> arma::mat {
  arma::rowvec v = h(X, theta_) % (1. - h(X, theta_)); // Variance function
  return X * arma::diagmat(v) * X.t();
}

auto &LogisticRegression::get_theta() -> arma::rowvec {
  return theta_;
}

auto LogisticRegression::Wald(arma::mat &X) -> arma::vec {
  arma::vec v = arma::inv_sympd(hessian(X)).eval().diag();
  arma::vec p(theta_.n_elem, arma::fill::zeros);

  boost::math::chi_squared chisq(1);

  for(arma::uword i = 0; i < theta_.n_elem; i++) {
    p[i] = boost::math::cdf(boost::math::complement(chisq, (theta_[i] * theta_[i]) / v[i]));
  }

  return p;
}
