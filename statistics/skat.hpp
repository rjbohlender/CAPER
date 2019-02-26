//
// Created by Bohlender,Ryan James on 8/22/18.
//

#ifndef PERMUTE_ASSOCIATE_SKAT_HPP
#define PERMUTE_ASSOCIATE_SKAT_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>

#include "../data/gene.hpp"
#include "../data/covariates.hpp"

// Unadjusted Structs
struct SKATParam {
  explicit SKATParam(arma::mat &Z1);

  // Fields
  double MuQ;
  double VarQ;
  double KerQ;
  arma::vec lambda;
  double VarRemain;
  double Df;
  arma::vec tau;

  static constexpr double rho_[8] = {0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 1};
};

struct LiuParam {
  LiuParam() = default;
  explicit LiuParam(arma::vec &c1);
  LiuParam(arma::vec &lambda, bool mod_lambda);

  // Fields
  double muQ;
  double sigmaQ;
  double varQ;
  double s1;
  double s2;
  double beta1;
  double beta2;

  bool type1;

  double muX;
  double sigmaX;

  double a, d, l;
};

struct PvalueLambda {
  PvalueLambda(arma::vec &lambda, double Q);

  double Get_Liu_Pval_MOD_Lambda(double Q, arma::vec &lambda, bool log_p = false);
  std::string Get_Liu_Pval_MOD_Lambda_Zero(double Q, LiuParam &param);

  // Fields
  double p_val;
  double p_val_liu;

  bool is_converged;

  double p_val_log;
  std::string p_val_zero_msg;
};

struct SKAT_Davies {
  SKAT_Davies(double Q, arma::vec &lambda);

  int ifault;
  double res;
};

struct SKAT_Optimal_GetQ {
  SKAT_Optimal_GetQ(arma::mat &Z,
					arma::vec &res,
					arma::vec &rall,
					arma::mat &res_out,
					int nResampling = 0);

  arma::vec Q_r;
  arma::mat Q_sim;
};

struct SKAT_EachQ {
  SKAT_EachQ(arma::vec &Qall, std::vector<arma::vec> &lambda_all);

  arma::vec pval;
  arma::vec pmin_q;
  double pmin;
};

struct SKAT_Integrate_Davies {
  SKAT_Integrate_Davies(arma::vec &pmin_q, SKATParam &param_m, arma::vec &rall);

  double operator()(double x);

  // fields
  const arma::uword nr = 8;
  arma::vec pmin_q;
  SKATParam param_m;
  arma::vec rall;

};

struct SKAT_Integrate_Liu {
  SKAT_Integrate_Liu(arma::vec &pmin_q, SKATParam &param_m, arma::vec &rall);

  double operator()(double x);

  // fields
  const arma::uword nr = 8;
  arma::vec pmin_q;
  SKATParam param_m;
  arma::vec rall;

};

double SKAT_Optimal_Pvalue_Davies(arma::vec &pmin_q, SKATParam &param_m, arma::vec &rall, double pmin);

#endif //PERMUTE_ASSOCIATE_SKAT_HPP
