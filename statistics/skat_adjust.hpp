//
// Created by Bohlender,Ryan James on 8/23/18.
//

#ifndef PERMUTE_ASSOCIATE_SKAT_ADJUST_HPP
#define PERMUTE_ASSOCIATE_SKAT_ADJUST_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <memory>

#include "skat.hpp"
#include "../data/gene.hpp"
#include "../data/covariates.hpp"

arma::vec Get_Lambda(arma::mat &K);

// Adjusted
struct SKAT_Residuals_Logistic {
  SKAT_Residuals_Logistic() = default;
  SKAT_Residuals_Logistic(Covariates &cov, int nresampling);
  SKAT_Residuals_Logistic(arma::mat &cov, arma::vec &mu, arma::vec &Y, int nresampling);

  void calcX1(arma::mat &cov);

  arma::uword nResampling = 0;
  arma::mat X1;
  arma::vec res;
  arma::mat res_out;
  arma::vec mu;
  arma::vec pi_1;
};

// Small sample adjustment
struct SKAT_Adjust {
  SKAT_Adjust(Gene &gene, Covariates &cov, const std::string &k, const std::string &kernel, int a, int b, std::shared_ptr<SKAT_Residuals_Logistic> &re2_, std::map<std::string, std::shared_ptr<SKAT_Optimal_GetQ>> &Q_sim_all);

  arma::vec weights;

  SKAT_Residuals_Logistic re1;

  double p_value = NAN;
};

// SKAT_Logistic_VarMatching_GetParam
struct SKAT_Logistic_VarMatching_Param {
  SKAT_Logistic_VarMatching_Param() = default;
  SKAT_Logistic_VarMatching_Param(arma::mat &Z1, arma::vec &p_all, arma::vec &Q_sim);
  SKAT_Logistic_VarMatching_Param(const SKAT_Logistic_VarMatching_Param &other) = default;
  SKAT_Logistic_VarMatching_Param(SKAT_Logistic_VarMatching_Param &&other) = default;

  SKAT_Logistic_VarMatching_Param &operator=(const SKAT_Logistic_VarMatching_Param &rhs) = default;

  void calc_param(arma::vec &lambda, arma::mat &U, arma::vec &p_all, arma::vec &Q_sim);

  double SKAT_Get_Var_Elements(arma::vec &m4,
                               arma::vec &p_all,
                               arma::subview_col<double> Ui,
                               arma::subview_col<double> Uj);
  double SKAT_Get_DF_Sim(arma::vec &Q_sim);
  double SKAT_Get_Kurtosis(arma::vec &x);
  double SKAT_Get_Skewness(arma::vec &x);
  void only_sim(arma::mat &Z1, arma::vec &p_all, arma::vec &Q_sim);

  std::pair<arma::vec, arma::mat>  Get_Lambda_U_From_Z(arma::mat &Z1);

  double muQ;
  double varQ;
  double df;
  arma::vec zeta;
  arma::vec var_i;
  arma::vec lambda_new;
  LiuParam param_noadj;
  arma::uword nlambda;
};

struct SKAT_Optimal_Logistic_VarMatching {
  SKAT_Optimal_Logistic_VarMatching() = default;
  SKAT_Optimal_Logistic_VarMatching(SKAT_Residuals_Logistic &re1,
									SKAT_Residuals_Logistic &re2,
									std::shared_ptr<SKAT_Optimal_GetQ> &Q_sim_all,
									arma::sp_mat &Z,
									const std::string &kernel,
									arma::vec &weights);

  arma::vec pval_each;
  arma::vec qval_each;
  double minp;
  double rho_est;
  double pval;

  static constexpr std::array<double, 8> rcorr {0.0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 0.999}; // Not one to prevent errors
};

struct SKAT_Optimal_Param_VarMatching {
  SKAT_Optimal_Param_VarMatching() = default;
  SKAT_Optimal_Param_VarMatching(arma::mat &Z1, arma::vec &rall, arma::vec &p_all, arma::vec &res, arma::mat &res_moments);

  SKAT_Optimal_Param_VarMatching(const SKAT_Optimal_Param_VarMatching &other) = default;
  SKAT_Optimal_Param_VarMatching &operator=(const SKAT_Optimal_Param_VarMatching &rhs);


  arma::vec Q_sim;
  double VarRemain;
  arma::vec tau;

  SKAT_Logistic_VarMatching_Param param;
};

// SKAT_PValue_Logistic_VarMatching
struct SKAT_PValue_Logistic_VarMatching {
  SKAT_PValue_Logistic_VarMatching() = default;
  // Z1 is from Z2_all
  SKAT_PValue_Logistic_VarMatching(double Q, arma::mat &Z1, arma::vec &p_all, arma::subview_col<double> Q_sim);

  double p_value;
  double p_value_noadj;

  SKAT_Logistic_VarMatching_Param param;
};

// SKAT_Optimal_Each_Q_VarMatching
struct SKAT_Optimal_Each_Q_VarMatching {
  SKAT_Optimal_Each_Q_VarMatching() = default;
  SKAT_Optimal_Each_Q_VarMatching(SKAT_Optimal_Param_VarMatching &param_m,
                                    arma::vec &Q_all,
                                    arma::vec &rall,
                                    std::vector<arma::mat> &Z2_all,
                                    arma::vec &p_all,
                                    arma::mat &Q_sim_all);
  double SKAT_Optimal_Kurtosis_Mixture(double df1, double df2, double v1, double a1, double a2);

  arma::vec pval;
  arma::vec pmin_q;
  double pmin;
};

// SKAT_Optimal_Get_Pvalue_VarMatching
struct SKAT_Optimal_Get_Pvalue_VarMatching {
  SKAT_Optimal_Get_Pvalue_VarMatching(arma::vec &Q_all,
									  arma::mat &Z1,
									  arma::vec &rall,
									  arma::vec &p_all,
									  arma::mat &Q_sim_all,
									  arma::vec &res,
									  arma::mat &res_moments);

  SKAT_Optimal_Param_VarMatching param_m;
  SKAT_Optimal_Each_Q_VarMatching each_info;

  arma::vec pval_each;
  double pval;
  double SKAT_Optimal_Pvalue_VarMatching(arma::vec pmin_q,
                                         double muQ,
                                         double varQ,
                                         double df,
                                         arma::vec tau,
                                         arma::vec &rall,
                                         double pmin = NAN);
};
#endif //PERMUTE_ASSOCIATE_SKAT_ADJUST_HPP
