//
// Created by Bohlender,Ryan James on 8/1/18.
//

#ifndef PERMUTE_ASSOCIATE_CARVATASK_HPP
#define PERMUTE_ASSOCIATE_CARVATASK_HPP

#define ARMA_DONT_USE_WRAPPER

#include <string>
#include <cmath>
#include <unordered_map>
#include <cassert>
#include <armadillo>
#include <set>

#include <boost/optional.hpp>

#include "../data/gene.hpp"
#include "../data/covariates.hpp"
#include "../statistics/methods.hpp"
#include "../data/result.hpp"
#include "../data/permutation.hpp"

enum class Stage {
  Stage1,
  Stage2,
  Done,
  Power
};

class CARVATask {
public:
  Stage stage;
  Gene gene;
  std::shared_ptr<Covariates> cov;
  Methods methods;
  size_t success_threshold;
  size_t npermutations;
  size_t offset;
  size_t termination;
  std::vector<std::vector<int8_t>> &permutations;
  const TaskParams tp;
  std::unordered_map<std::string, Result> results;
  std::unordered_map<std::string, Permute> permute;

  // Constructors
  CARVATask(Stage stage_,
			Gene gene_,
			const std::shared_ptr<Covariates> &cov_,
			TaskParams tp_,
			arma::uword succ_thresh_,
			arma::uword nperm_,
			arma::uword offset_,
			arma::uword termination_,
			std::vector<std::vector<int8_t>> &perm_);
  CARVATask(Stage stage_,
			Gene &gene_,
			std::shared_ptr<Covariates> cov_,
			TaskParams tp_,
			std::vector<std::vector<int8_t>> &perm_);

  // Free memory
  void cleanup();

  Covariates &get_cov();
  int max_permutations();
  std::vector<std::vector<int8_t>> &get_permutations();

  void calc_multitranscript_pvalues();
  bool has_multiple_transcripts();
};

#endif //PERMUTE_ASSOCIATE_CARVATASK_HPP
