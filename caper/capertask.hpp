//
// Created by Bohlender,Ryan James on 8/1/18.
//

#ifndef PERMUTE_ASSOCIATE_CAPERTASK_HPP
#define PERMUTE_ASSOCIATE_CAPERTASK_HPP

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif

#include <armadillo>
#include <cassert>
#include <cmath>
#include <set>
#include <string>
#include <unordered_map>

#include <boost/optional.hpp>

#include "../data/covariates.hpp"
#include "../data/gene.hpp"
#include "../data/permutation.hpp"
#include "../data/result.hpp"
#include "../statistics/methods.hpp"

enum class Stage { Stage1, Stage2, Done, Power };

class CAPERTask {
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
  CAPERTask(Stage stage_,
            Gene gene_,
            const std::shared_ptr<Covariates> &cov_,
            TaskParams tp_,
            arma::uword succ_thresh_,
            arma::uword nperm_,
            arma::uword offset_,
            arma::uword termination_,
            std::vector<std::vector<int8_t>> &perm_);
  CAPERTask(Stage stage_,
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

#endif // PERMUTE_ASSOCIATE_CAPERTASK_HPP
