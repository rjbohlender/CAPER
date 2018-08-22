//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_METHODS_HPP
#define PERMUTE_ASSOCIATE_METHODS_HPP

#include <map>
#include <string>
#include <armadillo>

#include "../data/gene.hpp"
#include "../data/covariates.hpp"

arma::vec rank(arma::vec &v, const char *direction);

enum class Kernel {
  Linear,
  wLinear,
  Quadratic,
  IBS,
  wIBS,
  twoWayX
};

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

  static constexpr double rho_[8] = {0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 1};;
};

struct LiuParam {
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
  SKAT_Optimal_GetQ(arma::mat &Z1, arma::vec &res, arma::vec &rall);

  arma::vec Q_r;
  arma::vec Q_sim;
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

class Methods {
public:
  explicit Methods(std::string method);
  Methods(std::string method, std::string kernel);

  double call(arma::mat &Xmat, arma::vec &Yvec);
  double call(arma::mat &Xmat, arma::vec &Yvec, arma::vec &weights);
  double call(arma::mat &Xmat, arma::vec &Yvec, arma::vec &weights, bool score_only_minor, double site_penalty);
  double call(const std::string &k, Gene &gene, Covariates &cov);
  double call(const std::string &k, Gene &gene, Covariates &cov, arma::colvec &weights);
  double call(const std::string &k, Gene &gene, Covariates &cov, bool shuffle);

  std::string str();

  void clear(std::vector<std::string> &v);

private:
  const std::string method_;

  // SKAT support fields
  Kernel kernel_;
  std::map<std::string, arma::mat> K_;
  // Weights for SKAT-O
  static constexpr double rho_[8] = {0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 1};

  double BURDEN(arma::mat &Xmat, Covariates &cov, arma::vec &weights);
  double CALPHA(arma::mat &Xmat, arma::vec &Yvec);
  double CMC(arma::mat &Xmat, arma::vec &Yvec, double maf = 0.05);
  double SKAT(arma::mat &Xmat,
			  Covariates &cov,
			  arma::vec &weights,
			  const std::string &k,
			  bool shuffle = false,
			  int a = 1,
			  int b = 25);
  double SKATO(arma::mat &Xmat,
			   Covariates &cov,
			   arma::vec &weights,
			   const std::string &k,
			   bool shuffle = false,
			   int a = 1,
			   int b = 25);
  double WSS(arma::mat &Xmat, arma::colvec &Yvec);
  double VAAST(arma::mat &Xmat,
			   arma::colvec &Yvec,
			   arma::colvec &log_casm,
			   bool score_only_minor = true,
			   bool score_only_alternative = true,
			   double site_penalty = 2.0);
  double VT(arma::mat &Xmat, arma::colvec &Yvec);

  // VAAST support member functions
  arma::vec LRT(arma::vec &case_allele1,
				arma::vec &control_allele1,
				arma::vec &case_allele0,
				arma::vec &control_allele0);
  arma::vec log_likelihood(arma::vec &freq, arma::vec &allele0, arma::vec &allele1);

  // Kernel member functions
  arma::mat kernel_Linear(arma::mat &Xmat);
  arma::mat kernel_wLinear(arma::mat &Xmat, arma::vec &weights);
  arma::mat kernel_IBS(arma::mat &Xmat, arma::uword &n, arma::uword &p);
  arma::mat kernel_wIBS(arma::mat &Xmat, arma::uword &n, arma::uword &p, arma::vec &weights);
  arma::mat kernel_Quadratic(arma::mat &Xmat);
  arma::mat kernel_twoWayX(arma::mat &Xmat, arma::uword n, arma::uword p);

  // SKAT-O support
  double SKAT_Optimal_Pvalue_Davies(arma::vec &pmin_q, SKATParam &param_m, arma::vec &rall, double pmin);
};

#endif //PERMUTE_ASSOCIATE_METHODS_HPP
