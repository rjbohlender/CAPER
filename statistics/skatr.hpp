//
// Created by Bohlender,Ryan James on 9/22/18.
//

#ifndef PERMUTE_ASSOCIATE_SKATR_HPP
#define PERMUTE_ASSOCIATE_SKATR_HPP

#include "../data/covariates.hpp"

class SKATR_Null {
public:
  explicit SKATR_Null(Covariates &cov);

  // Handle permutation
  auto shuffle(bool skip_svd = false) noexcept -> void;

  auto get_U0() noexcept -> arma::vec;
  auto get_pi0() noexcept -> arma::vec;
  auto get_Yv() noexcept -> arma::vec;
  auto get_Yh() noexcept -> arma::vec;
  auto get_Ux() noexcept -> const arma::mat &;
  auto get_coef() noexcept -> arma::rowvec;
  auto get_Y() noexcept -> arma::vec;
  auto get_X() noexcept -> const arma::mat &;

private:
  CRandomMersenne crand;
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

class SKATR_Linear_Null {
public:
  explicit SKATR_Linear_Null(Covariates &cov);

  // Handle permutation
  auto shuffle() noexcept -> void;

  auto get_U0() noexcept -> arma::vec;
  auto get_s2() noexcept -> double;
  auto get_Ux() noexcept -> const arma::mat &;
  auto get_coef() noexcept -> const arma::rowvec &;
  auto get_Y() noexcept -> arma::vec;
  auto get_X() noexcept -> const arma::mat &;

private:
  CRandomMersenne crand;
  arma::vec U0;          // Residuals
  double s2;             // The residual variance
  arma::mat Ux;          // ;
  arma::rowvec coef;     // model coefficients -- theta_
  arma::vec Y;           // Phenotype
  arma::mat X;           // Design matrix
  arma::uvec indices;    // Index vector
};

#endif //PERMUTE_ASSOCIATE_SKATR_HPP
