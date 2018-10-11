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

#include <boost/optional.hpp>

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

struct TaskParams {
  // Permutation paramters
  arma::uword success_threshold;

  arma::uword stage_1_permutations;
  arma::uword stage_2_permutations;
  arma::uword total_permutations;

  // For SKATO, SKAT, BURDEN
  bool alternate_permutation;

  // Method
  std::string method;

  // General options
  std::string program_path;
  std::string genotypes_path;
  std::string covariates_path;
  std::string ped_path;

  boost::optional<std::string> bed;
  boost::optional<std::string> weight;

  bool verbose;

  // Output permutations
  boost::optional<std::string> permute_set;

  // Detailed VAAST output
  std::string full_command;
  std::string output_path;

  // VAAST
  arma::uword group_size;
  bool score_only_minor;
  bool score_only_alternative;
  bool testable;

  // CMC
  double maf;

  size_t nthreads;

  // Gene list
  boost::optional<std::string> gene_list;

  // SKAT Parameters
  std::string kernel; // Kernel selection
  bool adjust; // Sample size adjustment
  double a; // Beta weight parameters
  double b; // Beta weight parameters
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
		   TaskParams &tp,
		   int succ_thresh,
		   int s1_perm,
		   int s2_perm,
		   std::vector<std::vector<int32_t>> &perm);
  TaskArgs(Stage stage,
		   Gene gene,
		   Covariates cov,
		   TaskParams &tp,
		   std::vector<std::vector<int32_t>> &perm);
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
  TaskParams &get_tp();

  // Setter
  void set_stage(Stage stage);

  void calc_multitranscript_pvalues();
  bool has_multiple_transcripts();

  Permute &get_permute(const std::string &k);

  int get_a();
  int get_b();

private:
  Stage stage_;
  Gene gene_;
  Covariates cov_;
  Methods method_;
  TaskParams tp_;

  int stage_1_permutations_;
  int stage_2_permutations_;
};

#endif //PERMUTE_ASSOCIATE_TASKARGS_HPP
