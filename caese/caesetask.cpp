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
  : gene(gene), cov(cov), methods(tp, cov), tp(tp), permutations_(tp.nperm) {
  for (const auto &k : gene.get_transcripts()) {
	results[k] = Result(gene.gene_name, k, !gene.is_polymorphic(k));
	if(tp.seed) {
	  permute_[k] = Permute(*tp.seed);
	} else {
	  permute_[k] = Permute();
	}
  }
}
CAESETask::CAESETask(Stage stage,
					 Gene gene,
					 const std::shared_ptr<Covariates> &cov,
					 TaskParams &tp,
					 std::vector<std::vector<int32_t>> &perm)
	: gene(gene), cov(cov), methods(tp, cov), tp(tp), permutations_(tp.nperm) {
  for (const auto &k : gene.get_transcripts()) {
	results[k] = Result(gene.gene_name, k, !gene.is_polymorphic(k));
	if(tp.seed) {
	  permute_[k] = Permute(*tp.seed);
	} else {
	  permute_[k] = Permute();
	}
  }
}

auto CAESETask::cleanup() -> void {

}

auto CAESETask::get_cov() -> Covariates & {
  return *cov;
}

auto CAESETask::get_methods() -> Methods & {
  return methods;
}

auto CAESETask::get_tp() -> TaskParams & {
  return tp;
}

auto CAESETask::get_permute(const std::string &k) -> Permute & {
  return permute_[k];
}

auto CAESETask::get_npermutations() -> int {
  return permutations_;
}
