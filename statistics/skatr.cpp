//
// Created by Bohlender,Ryan James on 9/22/18.
//

#include "skatr.hpp"

SKATR_Null::SKATR_Null(Covariates cov)
: crand(std::random_device{}()), cov_(cov) {
  X = cov_.get_covariate_matrix();
  Y = cov_.get_phenotype_vector();
  pi0 = cov_.get_fitted();
  Yv = pi0 % (1 - pi0);
  Yh = arma::sqrt(Yv);

  arma::mat U, V;
  arma::vec S;
  arma::svd_econ(U, S, V, arma::diagmat(Yh) * X);

  Ux = arma::diagmat(Yh) * U;
  U0 = Y - pi0;
}

auto SKATR_Null::shuffle(arma::vec &phenotypes) -> void {
  cov_.set_phenotype_vector(phenotypes);
  cov_.refit_permuted();

  pi0 = cov_.get_fitted();
  Yv = pi0 % (1 - pi0);
  Yh = arma::sqrt(Yv);
  Y = phenotypes;

  arma::mat U, V;
  arma::vec S;
  bool success = arma::svd_econ(U, S, V, arma::diagmat(Yh) * X, "left");
  if (!success) {
    std::cerr << "Yh: " << Yh.t();
    std::cerr << "Yv: " << Yv.t();
    std::cerr << "Y: " << Y.t();
    for (int i = 0; i < Y.n_elem; i++) {
      std::cerr << Y[i];
      for (int j = 0; j < X.n_cols; j++) {
        std::cerr << "\t" << X(i, j);
      }
      std::cerr << std::endl;
    }
  }

  Ux = arma::diagmat(Yh) * U;
  U0 = Y - pi0;
}

auto SKATR_Null::get_U0() noexcept -> arma::vec {
  return U0;
}

auto SKATR_Null::get_pi0() noexcept -> arma::vec {
  return pi0;
}

auto SKATR_Null::get_Yv() noexcept -> arma::vec {
  return Yv;
}

auto SKATR_Null::get_Yh() noexcept -> arma::vec {
  return Yh;
}

auto SKATR_Null::get_Ux() noexcept -> const arma::mat & {
  return Ux;
}

auto SKATR_Null::get_coef() noexcept -> arma::rowvec {
  return coef;
}

auto SKATR_Null::get_Y() noexcept -> arma::vec {
  return Y;
}

auto SKATR_Null::get_X() noexcept -> const arma::mat & {
  return X;
}

SKATR_Linear_Null::SKATR_Linear_Null(Covariates cov)
: cov_(cov) {
  X = cov_.get_covariate_matrix();
  Y = cov_.get_phenotype_vector();

  arma::mat U, V;
  arma::vec s;
  arma::svd_econ(U, s, V, X.t());

  Ux = U;
  U0 = cov_.get_residuals();
  s2 = arma::accu(arma::pow(U0, 2)) / (Y.n_elem - X.n_rows); // Residual sum of squares over residual df
}

auto SKATR_Linear_Null::shuffle(arma::vec &phenotypes) -> void {
  cov_.set_phenotype_vector(phenotypes);
  cov_.refit_permuted();

  Y = phenotypes;

  U0 = cov_.get_residuals();
  s2 = arma::accu(arma::pow(U0, 2)) / (Y.n_elem - X.n_rows);
}

auto SKATR_Linear_Null::get_U0() noexcept -> arma::vec {
  return U0;
}

auto SKATR_Linear_Null::get_s2() noexcept -> double {
  return s2;
}

auto SKATR_Linear_Null::get_Ux() noexcept -> const arma::mat & {
  return Ux;
}

auto SKATR_Linear_Null::get_coef() noexcept -> const arma::rowvec & {
  return coef;
}

auto SKATR_Linear_Null::get_Y() noexcept -> arma::vec {
  return Y;
}
auto SKATR_Linear_Null::get_X() noexcept -> const arma::mat & {
  return X;
}
