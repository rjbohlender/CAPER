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
  arma::rowvec &get_theta();

  arma::vec Wald(arma::mat &X);

private:
  arma::colvec odds_;
  arma::rowvec theta_; // Vector of weights - Beta on wikipedia
  arma::vec mu_; // Probability vector - fitted.values in R
  arma::vec eta_;
  double mean_;

  arma::rowvec h(arma::mat &Xmat, arma::rowvec &t);
  double cost(arma::mat &Xmat, arma::colvec &Yvec, arma::rowvec &t);
  void train(arma::mat &Xmat, arma::colvec &Yvec);

  arma::rowvec IRLS_SVDNEWTON(arma::mat &Xmat, arma::colvec &Yvec);

  arma::mat hess(arma::mat &X);

  void calculate_odds(arma::mat &Xmat);
  void calculate_mean(arma::mat &Xmat);
  void calculate_probability(arma::mat &Xmat);
  void calculate_eta();
};

#endif //PERMUTE_ASSOCIATE_LOGISTIC_REGRESSION_HPP
