//
// Created by Bohlender,Ryan James on 8/2/18.
//

#ifndef PERMUTE_ASSOCIATE_LOGISTIC_REGRESSION_HPP
#define PERMUTE_ASSOCIATE_LOGISTIC_REGRESSION_HPP 1

#include <armadillo>
#include <cmath>
#include <vector>

class LogisticRegression {
public:
  LogisticRegression(arma::mat &Xmat, arma::colvec &Yvec);

  auto get_odds() -> arma::vec&;
  auto get_mean() -> double;
  auto get_probability() -> arma::vec&;
  auto get_eta() -> arma::vec&;
  auto get_theta() -> arma::rowvec&;

  auto Wald(arma::mat &X) -> arma::vec;

private:
  arma::colvec odds_;
  arma::rowvec theta_; // Vector of weights - Beta on wikipedia
  arma::vec mu_; // Probability vector - fitted.values in R
  arma::vec eta_;
  double mean_;

  auto h(arma::mat &Xmat, arma::rowvec &t) -> arma::rowvec;
  auto cost(arma::mat &Xmat, arma::colvec &Yvec, arma::rowvec &t) -> double;

  // Algorithms for finding the optimum
  auto gradient_descent(arma::mat &Xmat, arma::colvec &Yvec) -> arma::rowvec;
  auto irls_svdnewton(arma::mat &Xmat, arma::colvec &Yvec) -> arma::rowvec;

  auto hessian(arma::mat &X) -> arma::mat;

  auto calculate_odds(arma::mat &Xmat) -> void;
  auto calculate_mean(arma::mat &Xmat) -> void;
  auto calculate_probability(arma::mat &Xmat) -> void;
  auto calculate_eta() -> void;
};

#endif //PERMUTE_ASSOCIATE_LOGISTIC_REGRESSION_HPP
