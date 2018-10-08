//
// Created by Bohlender,Ryan James on 8/2/18.
//

#include <cassert>
#include <iostream>
#include "logistic_regression.hpp"

/// Armadillo stores matrices in column major format. Iteration proceeds along columns.

LogisticRegression::LogisticRegression(arma::mat &Xmat, arma::colvec &Yvec)
: theta_(arma::rowvec(Xmat.n_rows, arma::fill::ones)){
#ifndef NDEBUG
  std::cerr << "Cost: " << cost(Xmat, Yvec) << "\n";
  std::cerr << "h: " << h(Xmat) << "\n";
#endif

  // train(Xmat, Yvec);

  IRLS_SVDNEWTON(Xmat, Yvec);

  calculate_probability(Xmat);
  calculate_odds(Xmat);
  calculate_mean(Xmat);
  calculate_eta();


#ifndef NDEBUG
  std::cerr << "odds: " << odds_.t();
  std::cerr << "mu_: " << mu_.t();
  std::cerr << "eta_: " << eta_.t();
  std::cerr << "theta_: " << theta_;
  std::cerr << "mean: " << mean_ << "\n";
#endif
}

arma::rowvec LogisticRegression::h(arma::mat &Xmat) {
  return 1. / (1. + arma::exp(-(theta_ * Xmat)));
}

double LogisticRegression::cost(arma::mat &Xmat, arma::colvec &Yvec) {
  double m = Yvec.n_rows;

  arma::mat ret = -1. / m * (arma::log(h(Xmat)) * Yvec + arma::log(1 - h(Xmat)) * (1 - Yvec));

  assert(ret.n_rows == 1);
  assert(ret.n_cols == 1);

  return ret(0, 0);
}

void LogisticRegression::train(arma::mat &Xmat, arma::colvec &Yvec) {
  int iterations = 0;
  int max_iter = 1000000;
  double alpha = 0.1; // Learning rate
  double epochs = 10;
  double decay = 0.5 / epochs;
  double tol = 1e-10;
  double m = Xmat.n_cols;
  arma::rowvec grad;
  do {

    // Vectorized update
    grad = alpha * (Xmat * (h(Xmat).t() - Yvec)).t() / m;
    theta_ -= grad;

    alpha *= 1. / (1. + decay * iterations); // Learning rate update

#if 0
    std::cerr << "iteration: " << iterations << " theta: " << theta_;
    std::cerr << "cost: " << cost(Xmat, Yvec) << "\n";
#endif

    iterations++;
  } while(iterations < max_iter && arma::any(grad > tol));
#ifndef NDEBUG
  std::cerr << "Logistic regression iterations: " << iterations << "\n";
#endif
}

void LogisticRegression::calculate_odds(arma::mat &Xmat) {
  if(arma::any(mu_ <= 0)) {
    std::cerr << "mu_ contains zeros\n";
  }
  odds_ = mu_ / (1. - mu_);
}

arma::colvec &LogisticRegression::get_odds() {
  return odds_;
}

void LogisticRegression::calculate_mean(arma::mat &Xmat) {
  mean_  = arma::mean(mu_);
}

double LogisticRegression::get_mean() {
  return mean_;
}

arma::vec &LogisticRegression::get_probability() {
  return mu_;
}

void LogisticRegression::calculate_probability(arma::mat &Xmat) {
  // Prevent numerical issues
  mu_ = arma::clamp(h(Xmat).t(), 0.00001, 0.99999);
}

arma::vec &LogisticRegression::get_eta() {
  return eta_;
}

void LogisticRegression::calculate_eta() {
  eta_ = arma::log(mu_ / (1. - mu_));
}

/**
 * @brief Performs Iteratively Reweighted Least Squares
 * @param Xmat Features
 * @param Yvec Response vector
 *
 * Algorithm from: https://bwlewis.github.io/GLM/#svdnewton
 */
void LogisticRegression::IRLS_SVDNEWTON(arma::mat &Xmat, arma::colvec &Yvec) {
  // Iteration params
  double update;
  double tol = 1e-8;
  int iter = 0;
  int max_iter = 25;

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
    return arma::exp(x) / arma::pow((arma::exp(x) + 1), 2);
  };

  do {
    arma::vec gv = g(t);
    arma::vec varg = gv % (1 - gv);
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

  theta_ = (V * (arma::diagmat(1. / S) * (U.t() * t))).t();
}

arma::rowvec &LogisticRegression::get_theta() {
  return theta_;
}

