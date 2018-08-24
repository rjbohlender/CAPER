//
// Created by Bohlender,Ryan James on 8/1/18.
//

#ifndef PERMUTE_ASSOCIATE_TASKARGS_HPP
#define PERMUTE_ASSOCIATE_TASKARGS_HPP

#include <string>
#include <cmath>
#include <unordered_map>
#include <cassert>
#include <armadillo>

#include "gene.hpp"
#include "covariates.hpp"
#include "../statistics/methods.hpp"
#include "result.hpp"
#include "permutation.hpp"

enum class Stage {
  Stage1,
  Stage2,
  Done
};

class TaskArgs {
public:
  std::unordered_map<std::string, Result> results;
  std::vector<std::vector<int32_t>> &permutations;
  std::unordered_map<std::string, Permute> permute;
  int success_threshold;
  int stop_check_threshold;
  bool adjust;

  // Constructors
  TaskArgs(Stage stage,
		   Gene gene,
		   Covariates cov,
		   int succ_thresh,
		   int s1_perm,
		   int s2_perm,
		   const std::string &method,
		   const std::string &kernel,
		   std::vector<std::vector<int32_t>> &perm,
		   bool adjust);
  TaskArgs(const TaskArgs &ta);
  TaskArgs(TaskArgs &&ta) noexcept;
  TaskArgs &operator=(const TaskArgs &rhs);

  // Free memory
  void cleanup();

  // Getters
  const Stage &get_stage() const;
  Gene &get_gene();
  Covariates &get_cov();
  Methods &get_methods();
  int get_max_permutations();
  int get_npermutations();
  std::vector<std::vector<int32_t>> &get_permutations();
  Result &max_original_statistic();
  Result &min_empirical_pvalue();
  Result &min_mgit_p();
  Result &min_transcript_permutations();
  Result &min_transcript_successes();
  double get_mgit_pvalue(const std::string &k);
  int get_remaining();

  // Setter
  void set_stage(Stage stage);

  void calc_multitranscript_pvalues();
  bool has_multiple_transcripts();

  Permute &get_permute(const std::string &k);

private:
  Stage stage_;
  Gene gene_;
  Covariates cov_;
  Methods method_;

  int stage_1_permutations_;
  int stage_2_permutations_;
};

#endif //PERMUTE_ASSOCIATE_TASKARGS_HPP
