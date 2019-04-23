//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_METHODS_HPP
#define PERMUTE_ASSOCIATE_METHODS_HPP

#define ARMA_DONT_USE_WRAPPER

#include <map>
#include <string>
#include <armadillo>
#include <memory>
#include <tuple>

#include "../data/gene.hpp"
#include "../data/covariates.hpp"
#include "skat.hpp"
#include "skat_adjust.hpp"
#include "skatr.hpp"
#include "vt.hpp"

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
  Methods(TaskParams &tp, Covariates &cov);

  std::string str();

  void clear(std::vector<std::string> &v);

  // Wu, Guan, and Pankow 2016
  double BURDEN(Gene &gene, const std::string &k, bool shuffle, int a, int b);
  double CALPHA(Gene &gene, arma::vec &Y, const std::string &k);
  // Li and Leal 2008
  double CMC(Gene &gene, arma::vec &Y, const std::string &k, double maf = 0.005);
  // Morris and Zeggini 2010
  double RVT1(Gene &gene, arma::vec &Y, arma::mat &design, const std::string &k, bool linear);
  double RVT2(Gene &gene, arma::vec &Y, arma::mat &design, const std::string &k, bool linear);
  // Wu, Guan, and Pankow 2016
  double SKATR(Gene &gene,
			   const std::string &k,
			   bool shuffle,
			   int a,
			   int b,
			   bool detail,
			   bool linear,
			   bool permute);
  // Wu, Guan, and Pankow 2016
  double SKATRO(Gene &gene, const std::string &k, bool shuffle, int a, int b, bool detail = false, bool linear = false);
  double Vaast(Gene &gene,
			   arma::vec &Y,
			   const std::string &k,
			   bool score_only_minor = true,
			   bool score_only_alternative = true,
			   double site_penalty = 2.0,
			   arma::uword group_threshold = 4,
			   bool detail = false,
			   bool biallelic = false);
  double VT(Gene &gene, const std::string &k, bool shuffle);
  // Madsen, Browning 2009
  double WSS(Gene &gene, arma::vec &Y, const std::string &k);

private:
  const std::string method_;

  // SKAT support fields
  Kernel kernel_;
  std::map<std::string, arma::mat> K_;
  // Weights for SKAT-O
  static constexpr double rho_[8] = {0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 1};

#if 0
  // Wu et al. 2011
  double SKAT(arma::mat &Xmat,
			  Covariates &cov,
			  arma::vec &weights,
			  const std::string &k,
			  bool shuffle = false,
			  int a = 1,
			  int b = 25);
#endif
#if 0
  // Lee et al. 2012
  double SKATO(Gene &gene,
			   Covariates &cov,
			   arma::vec &weights,
			   const std::string &k,
			   bool shuffle = false,
			   int a = 1,
			   int b = 25,
			   bool adjust = true);
#endif

  // Check weighting
  void check_weights(Gene &gene, const std::string &k, int a = 1, int b = 25, bool no_weight = false);

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
  std::shared_ptr<SKATR_Linear_Null> lin_obj_;

  // VT Helper
  std::shared_ptr<VT_Res> vt_obj_;

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
