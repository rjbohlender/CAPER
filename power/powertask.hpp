//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#ifndef PERMUTE_ASSOCIATE_POWERTASK_HPP
#define PERMUTE_ASSOCIATE_POWERTASK_HPP

#include <string>
#include <armadillo>
#include "../utility/taskparams.hpp"
#include "../data/gene.hpp"
#include "../data/covariates.hpp"
#include "../statistics/methods.hpp"
#include "../carva/carvatask.hpp"

struct PowerRes {
  std::string gene;
  std::string transcript;
  std::string method;
  arma::uword ncases;
  arma::uword ncontrols;
  double successes;
  double bootstraps;
  double ratio;
  double alpha;
};

struct PowerTask {
  PowerTask(Stage stage,
			Gene gene,
			const std::shared_ptr<Covariates> &cov,
			TaskParams &tp,
			std::vector<std::vector<int32_t>> &perm);
  PowerTask(Stage stage,
			Gene gene,
			const std::shared_ptr<Covariates> &cov,
			TaskParams &tp,
			arma::uword succ_thresh,
			arma::uword s1_perm,
			arma::uword s2_perm,
			std::vector<std::vector<int32_t>> &perm);
  TaskParams &tp;
  Gene gene;
  std::shared_ptr<Covariates> cov;
  arma::uword nreps;
  Methods method;

  // Power Results
  std::vector<PowerRes> power_res_;
};

#endif //PERMUTE_ASSOCIATE_POWERTASK_HPP
