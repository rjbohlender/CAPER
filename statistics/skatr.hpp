//
// Created by Bohlender,Ryan James on 9/22/18.
//

#ifndef PERMUTE_ASSOCIATE_SKATR_HPP
#define PERMUTE_ASSOCIATE_SKATR_HPP

#include "../data/covariates.hpp"

class SKATR_Null {
public:
  SKATR_Null(Covariates &cov);

  // Handle permutation
  void shuffle();

  // Getters
  arma::vec get_U0();
  arma::vec get_pi0();
  arma::vec get_Yv();
  arma::vec get_Yh();
  arma::mat get_Ux();
  arma::vec get_coef();
  arma::vec get_Y();
  arma::mat get_X();

private:
  arma::vec U0;          // Residuals
  arma::vec pi0;         // The fitted probabilities
  arma::vec Yv;          // The variance vector
  arma::vec Yh;          // sqrt of the variance vector
  arma::mat Ux;          // ;
  arma::rowvec coef;     // model coefficients -- theta_
  arma::vec Y;           // Case-control status
  arma::mat X;           // Design matrix
  arma::uvec indices;    // Index vector
};

#endif //PERMUTE_ASSOCIATE_SKATR_HPP
