//
// Created by Bohlender,Ryan James on 2019-03-19.
//

#ifndef PERMUTE_ASSOCIATE_POWER_HPP
#define PERMUTE_ASSOCIATE_POWER_HPP

#include <armadillo>
#include "../data/gene.hpp"
#include "../data/covariates.hpp"
#include "../data/permutation.hpp"
#include "methods.hpp"

struct PowerRes {
  std::string gene;
  std::string transcript;
  std::string method;
  arma::uword ncases;
  arma::uword ncontrols;
  double successes;
  double bootstraps;
  double ratio;
  double alpha;
};

class Power {
public:
  Power(Methods &method, Gene &gene, Covariates &cov, TaskParams &tp, arma::uword nreps);

  std::vector<PowerRes> get_results();
private:
  // Paramters
  TaskParams tp_;

  // PRNG
  std::mt19937 gen_;

  // Need to modify the data matrix so we maintain a copied version
  Gene bootstrapped_;
  Covariates cov_;
  arma::vec phenotypes_;
  Permute permute_;

  // Sample information
  arma::uvec cases_;
  arma::uvec controls_;

  arma::uvec case_idx_;
  arma::uvec control_idx_;

  // Power Results
  std::vector<PowerRes> power_res_;

  // Result storage
  std::map<std::string, std::vector<arma::vec>> success_map_;

  // Sample data with replacement
  arma::sp_mat sample(arma::sp_mat &X, arma::uword ncases, arma::uword ncontrols);

  // Call method for power analysis
  double call_method(Methods &method, Gene &gene, Covariates &cov, arma::vec &phenotypes, TaskParams &tp, const std::string &k);
};

#endif //PERMUTE_ASSOCIATE_POWER_HPP
