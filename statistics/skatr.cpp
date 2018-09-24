//
// Created by Bohlender,Ryan James on 9/22/18.
//

#include "skatr.hpp"

SKATR_Null::SKATR_Null(Covariates &cov) {
  X = cov.get_covariate_matrix();
  Y = cov.get_phenotype_vector();
  pi0 = cov.get_probability();
  Yv = pi0 % (1 - pi0);
  Yh = arma::sqrt(Yv);

  arma::mat U, V;
  arma::vec S;
  arma::svd_econ(U, S, V, arma::diagmat(Yh) * X.t());

  Ux = arma::diagmat(Yh) * U;
  U0 = Y - pi0;

  indices = arma::regspace<arma::uvec>(0, Y.n_rows - 1);
}

void SKATR_Null::shuffle() {
  indices = arma::shuffle(indices);
}

arma::vec SKATR_Null::get_U0() {
  return U0;
}

arma::vec SKATR_Null::get_pi0() {
  return pi0(indices);
}

arma::vec SKATR_Null::get_Yv() {
  return Yv(indices);
}

arma::vec SKATR_Null::get_Yh() {
  return Yh(indices);
}

arma::mat SKATR_Null::get_Ux() {
  return Ux;
}

arma::vec SKATR_Null::get_coef() {
  return coef;
}

arma::vec SKATR_Null::get_Y() {
  return Y(indices);
}

arma::mat SKATR_Null::get_X() {
  return X;
}
