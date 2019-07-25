//
// Created by Bohlender,Ryan James on 10/31/18.
//

#include "vt.hpp"

#include <ctime>

VT_Res::VT_Res(const Covariates &cov)
: indices_(arma::regspace<arma::uvec>(0, cov.get_nsamples() - 1))
, cov_(cov) {}

auto VT_Res::shuffle(arma::vec &phenotypes) -> void {
  cov_.set_phenotype_vector(phenotypes);
  cov_.refit_permuted();
}

auto VT_Res::get_residuals() -> arma::vec {
  return cov_.get_residuals();
}
