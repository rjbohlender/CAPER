//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#include "powertask.hpp"

PowerTask::PowerTask(Stage stage, Gene gene, std::shared_ptr<Covariates> cov,
                     TaskParams &tp,
                     std::vector<std::vector<int8_t>> &perm)
    : tp(tp), gene(gene), cov(cov), nreps(tp.bootstrap_reps),
      method(tp, cov) {}

PowerTask::PowerTask(Stage stage, Gene gene,
                     const std::shared_ptr<Covariates> &cov, TaskParams &tp,
                     arma::uword succ_thresh, arma::uword nperm_,
                     arma::uword offset_, arma::uword termination_,
                     std::vector<std::vector<int8_t>> &perm)
    : tp(tp), gene(gene), cov(cov), nreps(tp.bootstrap_reps),
      method(tp, cov) {}
