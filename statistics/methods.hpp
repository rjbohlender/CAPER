//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_METHODS_HPP
#define PERMUTE_ASSOCIATE_METHODS_HPP

#include <map>
#include <string>
#include <armadillo>
#include <memory>

#include "../data/gene.hpp"
#include "../data/covariates.hpp"
#include "skat.hpp"
#include "skat_adjust.hpp"
#include "skatr.hpp"

arma::vec rank(arma::vec &v, const char *direction);

enum class Kernel {
  Linear,
  wLinear,
  Quadratic,
  IBS,
  wIBS,
  twoWayX
};

class Methods {
public:
  explicit Methods(std::string method);
  Methods(std::string method, std::string kernel, Covariates &cov);

  double call(const std::string &k, Gene &gene, Covariates &cov);
  double call(const std::string &k, Gene &gene, Covariates &cov, bool shuffle, int a, int b);

  std::string str();

  void clear(std::vector<std::string> &v);

private:
  const std::string method_;

  // SKAT support fields
  Kernel kernel_;
  std::map<std::string, arma::mat> K_;
  // Weights for SKAT-O
  static constexpr double rho_[8] = {0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 1};

  // Wu, Guan, and Pankow 2016
  double BURDEN(Gene &gene, SKATR_Null &obj, const std::string &k, bool shuffle, int a, int b);
  double CALPHA(arma::mat &Xmat, arma::vec &Yvec);
  double CMC(arma::mat &Xmat, arma::vec &Yvec, double maf = 0.05);
  // Wu et al. 2011
  double SKAT(arma::mat &Xmat,
			  Covariates &cov,
			  arma::vec &weights,
			  const std::string &k,
			  bool shuffle = false,
			  int a = 1,
			  int b = 25);
  // Wu, Guan, and Pankow 2016
  double SKATR(Gene &gene, SKATR_Null &obj, const std::string &k, bool shuffle, int a, int b);
  // Lee et al. 2012
  double SKATO(Gene &gene,
			   Covariates &cov,
			   arma::vec &weights,
			   const std::string &k,
			   bool shuffle = false,
			   int a = 1,
			   int b = 25,
			   bool adjust = true);
  // Wu, Guan, and Pankow 2016
  double SKATRO(Gene &gene, SKATR_Null &obj, const std::string &k, bool shuffle, int a=1, int b=25);
  double WSS(arma::mat &Xmat, arma::colvec &Yvec);
  double VAAST(Gene &gene,
				 Covariates &cov,
				 const std::string &k,
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

  // Check weighting
  void check_weights(Gene &gene, const std::string &k, int a=1, int b=25);

  // Kernel member functions
  arma::mat kernel_Linear(arma::mat &Xmat);
  arma::mat kernel_wLinear(arma::mat &Xmat, arma::vec &weights);
  arma::mat kernel_IBS(arma::mat &Xmat, arma::uword &n, arma::uword &p);
  arma::mat kernel_wIBS(arma::mat &Xmat, arma::uword &n, arma::uword &p, arma::vec &weights);
  arma::mat kernel_Quadratic(arma::mat &Xmat);
  arma::mat kernel_twoWayX(arma::mat &Xmat, arma::uword n, arma::uword p);

  // SKATO Support
  std::shared_ptr<SKAT_Residuals_Logistic> re2;
  std::map<std::string, std::shared_ptr<SKAT_Optimal_GetQ>> Q_sim_all;

  // SKATR Null Model
  std::shared_ptr<SKATR_Null> obj_;

  // Davies method
  double SKAT_pval(double Q, arma::vec lambda);

  // Liu Method
  double Liu_qval_mod(double pval, arma::vec lambda);
  double Liu_pval(double Q, arma::vec lambda);
  double Saddlepoint(double Q, arma::vec lambda);

  template<class T>
  int sgn(T x);
};

#endif //PERMUTE_ASSOCIATE_METHODS_HPP
