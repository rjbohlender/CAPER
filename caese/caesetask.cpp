//
// Created by Bohlender,Ryan James on 2019-06-06.
//

#include "caesetask.hpp"

CAESETask::CAESETask(Stage stage,
					 Gene gene,
					 const std::shared_ptr<Covariates> &cov,
					 TaskParams &tp,
					 arma::uword succ_thresh,
					 arma::uword s1_perm,
					 arma::uword s2_perm,
					 std::vector<std::vector<int32_t>> &perm)
  : gene_(gene), cov_(cov), methods_(tp, *cov), tp_(tp), permutations_(tp.total_permutations) {
  for (const auto &k : gene_.get_transcripts()) {
	results[k] = Result(gene_.get_gene(), k);
	permute_[k] = Permute();
  }
}
CAESETask::CAESETask(Stage stage,
					 Gene gene,
					 const std::shared_ptr<Covariates> &cov,
					 TaskParams &tp,
					 std::vector<std::vector<int32_t>> &perm)
	: gene_(gene), cov_(cov), methods_(tp, *cov), tp_(tp), permutations_(tp.total_permutations) {
  for (const auto &k : gene_.get_transcripts()) {
	results[k] = Result(gene_.get_gene(), k);
	permute_[k] = Permute();
  }
}

auto CAESETask::cleanup() -> void {

}

auto CAESETask::get_gene() -> Gene & {
  return gene_;
}

auto CAESETask::get_cov() -> Covariates & {
  return *cov_;
}

auto CAESETask::get_methods() -> Methods & {
  return methods_;
}

auto CAESETask::get_tp() -> TaskParams & {
  return tp_;
}

auto CAESETask::get_permute(const std::string &k) -> Permute & {
  return permute_[k];
}

auto CAESETask::get_npermutations() -> int {
  return permutations_;
}
