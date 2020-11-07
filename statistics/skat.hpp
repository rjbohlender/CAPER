//
// Created by Bohlender,Ryan James on 9/22/18.
//

#ifndef PERMUTE_ASSOCIATE_SKATR_HPP
#define PERMUTE_ASSOCIATE_SKATR_HPP

#include "../data/covariates.hpp"

// Davies method
double SKAT_pval(double Q, const arma::vec& lambda);

// Liu Method
double Liu_qval_mod(double pval, const arma::vec& lambda);
double Liu_pval(double Q, const arma::vec& lambda);
double Saddlepoint(double Q, const arma::vec& lambda);

template<class T>
int sgn(T x) {
  return (T(0) < x) - (x < T(0));
}

class SKATR_Null {
public:
  explicit SKATR_Null(Covariates cov);

  // Handle permutation
  auto shuffle(arma::vec &phenotypes) -> void;

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
  Covariates cov_;
};

class SKATR_Linear_Null {
public:
  explicit SKATR_Linear_Null(Covariates cov);

  // Handle permutation
  auto shuffle(arma::vec &phenotypes) -> void;

  auto get_U0() noexcept -> arma::vec;
  auto get_s2() noexcept -> double;
  auto get_Ux() noexcept -> const arma::mat &;
  auto get_coef() noexcept -> const arma::rowvec &;
  auto get_Y() noexcept -> arma::vec;
  auto get_X() noexcept -> const arma::mat &;

private:
  arma::vec U0;          // Residuals
  double s2;             // The residual variance
  arma::mat Ux;          // ;
  arma::rowvec coef;     // model coefficients -- theta_
  arma::vec Y;           // Phenotype
  arma::mat X;           // Design matrix
  Covariates cov_;
};

#endif //PERMUTE_ASSOCIATE_SKATR_HPP
