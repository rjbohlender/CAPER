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
  std::unordered_map<std::string, Result> results;
  std::vector<std::vector<int32_t>> &permutations;
  std::unordered_map<std::string, Permute> permute;
  int success_threshold;
  int stop_check_threshold;
  bool adjust;

  // Constructors
  CARVATask(Stage stage,
		   Gene gene,
		   const std::shared_ptr<Covariates> &cov,
		   TaskParams &tp,
		   arma::uword succ_thresh,
		   arma::uword s1_perm,
		   arma::uword s2_perm,
		   std::vector<std::vector<int32_t>> &perm);
  CARVATask(Stage stage,
		   Gene gene,
		   const std::shared_ptr<Covariates> &cov,
		   TaskParams &tp,
		   std::vector<std::vector<int32_t>> &perm);
  CARVATask(const CARVATask &ta);
  CARVATask(CARVATask &&ta) noexcept;
  CARVATask &operator=(const CARVATask &rhs);

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

private:
  Stage stage_;
  Gene gene_;
  std::shared_ptr<Covariates> cov_;
  Methods method_;
  TaskParams tp_;

  static const std::set<std::string> pvalue_methods_;

  int stage_1_permutations_;
  int stage_2_permutations_;
};

#endif //PERMUTE_ASSOCIATE_CARVATASK_HPP
