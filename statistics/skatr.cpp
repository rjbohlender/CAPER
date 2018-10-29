//
// Created by Bohlender,Ryan James on 9/22/18.
//

#include "skatr.hpp"

SKATR_Null::SKATR_Null(Covariates &cov)
: crand((int)time(nullptr)) {
  X = cov.get_covariate_matrix();
  Y = cov.get_phenotype_vector();
  pi0 = cov.get_fitted();
  Yv = pi0 % (1 - pi0);
  Yh = arma::sqrt(Yv);

  arma::mat U, V;
  arma::vec S;
  arma::svd_econ(U, S, V, arma::diagmat(Yh) * X.t());

  Ux = arma::diagmat(Yh) * U;
  U0 = Y - pi0;

  indices = arma::regspace<arma::uvec>(0, Y.n_rows - 1);

  std::cerr << "pi0: " << pi0.t();
  std::cerr << "Ux: " << Ux.row(0);
  std::cerr << "U0: " << U0.t();
  std::cerr << "Yv: " << Yv.t();
  std::cerr << "Yh: " << Yh.t();
}

auto SKATR_Null::shuffle() noexcept -> void{
  // Fisher-Yates shuffle
  for(arma::sword i = indices.n_elem - 1; i >= 0; i--) {
	auto j = static_cast<arma::sword>(crand.IRandom(0, i));
	// Swap
	arma::uword tmp = indices[j];
	indices[j] = indices[i];
	indices[i] = tmp;
  }

  arma::mat U, V;
  arma::vec s;
  arma::svd_econ(U, s, V, arma::diagmat(Yh(indices)) * X.t());

  Ux = arma::diagmat(Yh) * U;
}

auto SKATR_Null::get_U0() noexcept -> arma::vec {
  return U0(indices);
}

auto SKATR_Null::get_pi0() noexcept -> arma::vec {
  return pi0(indices);
}

auto SKATR_Null::get_Yv() noexcept -> arma::vec {
  return Yv(indices);
}

auto SKATR_Null::get_Yh() noexcept -> arma::vec {
  return Yh(indices);
}

auto SKATR_Null::get_Ux() noexcept -> const arma::mat & {
  return Ux;
}

auto SKATR_Null::get_coef() noexcept -> arma::rowvec {
  return coef;
}

auto SKATR_Null::get_Y() noexcept -> arma::vec {
  return Y(indices);
}

auto SKATR_Null::get_X() noexcept -> const arma::mat & {
  return X;
}

SKATR_Linear_Null::SKATR_Linear_Null(Covariates &cov)
: crand((int)time(nullptr)) {
  X = cov.get_covariate_matrix();
  Y = cov.get_phenotype_vector();

  arma::mat U, V;
  arma::vec s;
  arma::svd_econ(U, s, V, X.t());

  Ux = U;
  U0 = cov.get_residuals();
  s2 = arma::accu(arma::pow(U0, 2)) / (Y.n_elem - X.n_rows); // Residual sum of squares over residual df

  indices = arma::regspace<arma::uvec>(0, Y.n_rows - 1);
}

auto SKATR_Linear_Null::shuffle() noexcept -> void {
  // Fisher-Yates shuffle
  for(arma::sword i = indices.n_elem - 1; i >= 0; i--) {
    auto j = static_cast<arma::sword>(crand.IRandom(0, i));
    // Swap
    arma::uword tmp = indices[j];
    indices[j] = indices[i];
    indices[i] = tmp;
  }
  indices = arma::shuffle(indices);
}

auto SKATR_Linear_Null::get_U0() noexcept -> arma::vec {
  return U0(indices);
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
  return Y(indices);
}
auto SKATR_Linear_Null::get_X() noexcept -> const arma::mat & {
  return X;
}
