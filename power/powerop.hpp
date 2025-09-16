//
// Created by Bohlender,Ryan James on 2019-03-19.
//

#ifndef PERMUTE_ASSOCIATE_POWEROP_HPP
#define PERMUTE_ASSOCIATE_POWEROP_HPP

#include <armadillo>
#include "../data/gene.hpp"
#include "../data/covariates.hpp"
#include "../data/permutation.hpp"
#include "../statistics/methods.hpp"
#include "powertask.hpp"
#include "../utility/reporter.hpp"

class PowerOp {
public:
  PowerOp(PowerTask &pt, std::shared_ptr<PowerReporter> reporter, double seed, bool verbose);

  auto run() -> void;
  auto finish() -> void;
  auto is_done() const -> bool;
  auto get_task() -> PowerTask;

  bool done_;
  PowerTask caperTask;
private:
  // Paramters
  std::shared_ptr<PowerReporter> reporter_;

  auto power() -> void;

  // PRNG
  std::mt19937 gen_;

  // Need to modify the data matrix so we maintain a copied version
  Gene bootstrapped_;
  arma::vec phenotypes_;
  Permute permute_;

  // Sample information
  arma::uvec cases_;
  arma::uvec controls_;

  arma::uvec case_idx_;
  arma::uvec control_idx_;


  // Result storage
  std::map<std::string, std::vector<arma::vec>> success_map_;

  // Sample data with replacement
  arma::sp_mat sample(arma::sp_mat &X, arma::uword ncases, arma::uword ncontrols);

  // Call method for power analysis
  double call_method(Methods &method, Gene &gene, Covariates &cov, arma::vec &phenotypes, TaskParams &tp, const std::string &k);
};

#endif //PERMUTE_ASSOCIATE_POWEROP_HPP
