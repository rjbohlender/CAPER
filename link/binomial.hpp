//
// Created by Bohlender,Ryan James on 10/15/18.
//

#ifndef PERMUTE_ASSOCIATE_BINOMIAL_HPP
#define PERMUTE_ASSOCIATE_BINOMIAL_HPP

#include "link.hpp"

struct Binomial : Link {
  enum class LinkID {
    Logit,
    Probit,
    cloglog,
    Cauchit,
    Log
  };

  static const std::vector<std::string> links;
  const LinkID linkid;

  explicit Binomial(const std::string &link="logit");
  arma::vec link(arma::mat &X, arma::vec &beta) noexcept override;
  arma::vec link(arma::vec &mu) noexcept override;
  arma::vec linkinv(arma::mat &X, arma::vec &beta) noexcept override;
  arma::vec linkinv(arma::vec &eta) noexcept override;
  arma::vec variance(arma::vec &mu) noexcept override;
  arma::vec mueta(arma::vec &eta) noexcept override;

  LinkID check_linkid(const std::string &link);
};

#endif //PERMUTE_ASSOCIATE_BINOMIAL_HPP
