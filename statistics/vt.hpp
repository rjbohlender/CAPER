//
// Created by Bohlender,Ryan James on 10/31/18.
//

#ifndef PERMUTE_ASSOCIATE_VT_HPP
#define PERMUTE_ASSOCIATE_VT_HPP

#include "../data/covariates.hpp"

class VT_Res {
public:
  explicit VT_Res(const Covariates &cov);

  auto shuffle(arma::vec &phenotypes) -> void;

  auto get_residuals() -> arma::vec;

private:
  arma::uvec indices_;
  Covariates cov_;
};

#endif //PERMUTE_ASSOCIATE_VT_HPP
