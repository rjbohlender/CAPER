//
// Created by Bohlender,Ryan James on 10/31/18.
//

#ifndef PERMUTE_ASSOCIATE_VT_HPP
#define PERMUTE_ASSOCIATE_VT_HPP

#include "../data/covariates.hpp"
#include "../data/gene.hpp"

class VT_Res {
public:
  explicit VT_Res();
  auto initialize(Gene &gene, arma::vec &pheno, const std::string &ts) -> void;

  auto is_initialized(const std::string &ts) const -> bool;
  auto get_subset(const std::string &ts) -> arma::uvec&;
  auto get_geno(const std::string &ts) -> arma::vec&;
  auto sum_groups(arma::vec &X, arma::uvec &range, const std::string &ts) -> arma::vec;
  auto get_mCount(const std::string &ts) -> arma::vec&;
  auto get_oneToLen(const std::string &ts) -> arma::uvec&;
  auto get_csCountMeanpheno(const std::string &ts) -> arma::vec&;
  auto get_sqrtCsCountSquare(const std::string &ts) -> arma::vec&;

private:
  std::map<std::string, bool> initialized_;

  double meanpheno_;
  std::map<std::string, arma::vec> geno_;
  std::map<std::string, arma::uvec> snpid_;
  std::map<std::string, arma::uvec> subset_;
  std::map<std::string, arma::vec> count_;
  std::map<std::string, arma::vec> mCount_;
  std::map<std::string, arma::vec> countg_;
  std::map<std::string, arma::vec> countSquare_;
  std::map<std::string, arma::vec> weight_;
  std::map<std::string, arma::vec> countWeight_;
  std::map<std::string, arma::vec> group_;
  std::map<std::string, arma::vec> ugroup_;
  std::map<std::string, arma::vec> arr_;
  std::map<std::string, arma::uvec> oneToLen_;
  std::map<std::string, std::vector<arma::uvec>> whiches_;
  std::map<std::string, arma::vec> csCount_;
  std::map<std::string, arma::vec> csCountSquare_;
  std::map<std::string, arma::vec> csCountMeanpheno_;
  std::map<std::string, arma::vec> sqrtCsCountSquare_;
};

#endif //PERMUTE_ASSOCIATE_VT_HPP
