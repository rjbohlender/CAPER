//
// Created by Bohlender,Ryan James on 10/31/18.
//

#ifndef PERMUTE_ASSOCIATE_VT_HPP
#define PERMUTE_ASSOCIATE_VT_HPP

#include <armadillo>
#include <map>
#include <string>

#include "../data/gene.hpp"

class VT_Res {
public:
  explicit VT_Res();
  auto initialize(Gene &gene, const arma::vec &pheno, const std::string &ts)
      -> void;

  auto is_initialized(const std::string &ts) const -> bool;
  auto accumulate_thresholds(const std::string &ts, const arma::vec &values) const
      -> arma::vec;
  auto get_csCountMeanpheno(const std::string &ts) -> arma::vec &;
  auto get_sqrtCsCountSquare(const std::string &ts) -> arma::vec &;

private:
  std::map<std::string, bool> initialized_;
  std::map<std::string, arma::uvec> sorted_indices_;
  std::map<std::string, arma::uvec> threshold_offsets_;
  std::map<std::string, arma::vec> csCountMeanpheno_;
  std::map<std::string, arma::vec> sqrtCsCountSquare_;
};

#endif //PERMUTE_ASSOCIATE_VT_HPP
