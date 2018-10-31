//
// Created by Bohlender,Ryan James on 10/31/18.
//

#include "vt.hpp"

#include <ctime>

VT_Res::VT_Res(const Covariates &cov)
: residuals_(cov.get_residuals())
, indices_(arma::regspace<arma::uvec>(0, cov.get_nsamples() - 1))
, crand_((int)time(nullptr)) {}

auto VT_Res::shuffle() -> void {
  // Fisher-Yates shuffle
  for(arma::sword i = indices_.n_elem - 1; i >= 0; i--) {
	auto j = static_cast<arma::sword>(crand_.IRandom(0, i));
	// Swap
	arma::uword tmp = indices_[j];
	indices_[j] = indices_[i];
	indices_[i] = tmp;
  }
}

auto VT_Res::get_residuals() -> arma::vec {
  return residuals_(indices_);
}
