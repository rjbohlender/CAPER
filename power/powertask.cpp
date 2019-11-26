//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#include "powertask.hpp"

PowerTask::PowerTask(Stage stage,
					 Gene gene,
					 const std::shared_ptr<Covariates> &cov,
					 TaskParams &tp,
					 std::vector<std::vector<int32_t>> &perm)
	: tp(tp), gene(gene), cov(cov), method(tp, cov) {

}

PowerTask::PowerTask(Stage stage,
					 Gene gene,
					 const std::shared_ptr<Covariates> &cov,
					 TaskParams &tp,
					 arma::uword succ_thresh,
					 arma::uword s1_perm,
					 arma::uword s2_perm,
					 std::vector<std::vector<int32_t>> &perm)
	: tp(tp), gene(gene), cov(cov), method(tp, cov) {

}
