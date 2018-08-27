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

  train(Xmat, Yvec);
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
  return 1 / (1 + arma::exp(-(theta_ * Xmat)));
}

double LogisticRegression::cost(arma::mat &Xmat, arma::colvec &Yvec) {
  int m = Yvec.n_rows;

  arma::mat ret = -1. / m * (arma::log(h(Xmat)) * Yvec + arma::log(1 - h(Xmat)) * (1 - Yvec));

  assert(ret.n_rows == 1);
  assert(ret.n_cols == 1);

  return ret(0, 0);
}

void LogisticRegression::train(arma::mat &Xmat, arma::colvec &Yvec) {
  int iterations = 1000000;
  double alpha = 0.0005; // Learning rate
  double tol = 1e-9;
  double m = Xmat.n_cols;
  arma::rowvec grad;
  do {

    // Vectorized update
    grad = alpha * (Xmat * (h(Xmat).t() - Yvec)).t() / m;
    theta_ -= grad;

#if 0
    std::cerr << "iteration: " << iterations << " theta: " << theta_;
    std::cerr << "cost: " << cost(Xmat, Yvec) << "\n";
#endif

    iterations--;
  } while(iterations > 0 && arma::any(grad > tol));
#ifndef NDEBUG
  std::cerr << "Logistic regression iterations: " << iterations << "\n";
#endif
}

void LogisticRegression::calculate_odds(arma::mat &Xmat) {
  odds_ = mu_ / (1 - mu_);
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
  mu_ = h(Xmat).t();
}

arma::vec &LogisticRegression::get_eta() {
  return eta_;
}

void LogisticRegression::calculate_eta() {
  eta_ = arma::log(mu_ / (1. - mu_));
}

