//
// Created by Bohlender,Ryan James on 2019-06-06.
//

#ifndef PERMUTE_ASSOCIATE_CAESETASK_HPP
#define PERMUTE_ASSOCIATE_CAESETASK_HPP

#include <map>
#include <vector>

#include "../data/gene.hpp"
#include "../data/covariates.hpp"
#include "../data/permutation.hpp"
#include "../data/result.hpp"
#include "../carva/carvatask.hpp"

class CAESETask {
public:
  std::map<std::string, Result> results;
  // Constructors
  CAESETask(Stage stage,
			Gene gene,
			const std::shared_ptr<Covariates> &cov,
			TaskParams &tp,
			arma::uword succ_thresh,
			arma::uword s1_perm,
			arma::uword s2_perm,
			std::vector<std::vector<int32_t>> &perm);
  CAESETask(Stage stage,
			Gene gene,
			const std::shared_ptr<Covariates> &cov,
			TaskParams &tp,
			std::vector<std::vector<int32_t>> &perm);

  // Free memory
  auto cleanup() -> void;

  // Getters
  auto get_gene() -> Gene &;
  auto get_cov() -> Covariates &;
  auto get_methods() -> Methods &;
  auto get_tp() -> TaskParams &;
  auto get_npermutations() -> int;
  auto get_permute(const std::string &k) -> Permute &;

private:
  Gene gene_;
  std::shared_ptr<Covariates> cov_;
  Methods methods_;
  TaskParams tp_;

  const int permutations_;

  std::map<std::string, Permute> permute_;

};

#endif //PERMUTE_ASSOCIATE_CAESETASK_HPP
