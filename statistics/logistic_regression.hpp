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

  arma::colvec &get_odds();
  double get_mean();
  arma::vec &get_probability();
  arma::vec &get_eta();

private:
  arma::colvec odds_;
  arma::rowvec theta_; // Vector of weights - Beta on wikipedia
  arma::vec mu_; // Probability vector - fitted.values in R
  arma::vec eta_;
  double mean_;

  arma::rowvec h(arma::mat &Xmat);
  double cost(arma::mat &Xmat, arma::colvec &Yvec);
  void train(arma::mat &Xmat, arma::colvec &Yvec);

  void IRLS_SVDNEWTON(arma::mat &Xmat, arma::colvec &Yvec);

  void calculate_odds(arma::mat &Xmat);
  void calculate_mean(arma::mat &Xmat);
  void calculate_probability(arma::mat &Xmat);
  void calculate_eta();
};

#endif //PERMUTE_ASSOCIATE_LOGISTIC_REGRESSION_HPP
