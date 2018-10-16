//
// Created by Bohlender,Ryan James on 10/15/18.
//

#ifndef PERMUTE_ASSOCIATE_LINK_HPP
#define PERMUTE_ASSOCIATE_LINK_HPP

#include <armadillo>

struct Link {
  virtual arma::vec link(arma::mat &X, arma::vec &beta) noexcept = 0;
  virtual arma::vec link(arma::vec &mu) noexcept = 0;
  virtual arma::vec linkinv(arma::mat &X, arma::vec &beta) noexcept = 0;
  virtual arma::vec linkinv(arma::vec &eta) noexcept = 0;
  virtual arma::vec variance(arma::vec &mu) noexcept = 0;
  virtual arma::vec mueta(arma::vec &eta) noexcept = 0;
};

#endif //PERMUTE_ASSOCIATE_LINK_HPP
