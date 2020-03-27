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
#include "vt.hpp"

arma::vec rank(arma::vec &v, const char *direction);

enum class Kernel {
  Linear,
  wLinear
};

class Methods {
public:
  explicit Methods(std::string method);
  Methods(TaskParams &tp, const std::shared_ptr<Covariates> &cov);

  std::string str();

  void clear(std::vector<std::string> &v);

  // Wu, Guan, and Pankow 2016
  double BURDEN(Gene &gene, const std::string &k, arma::vec &phenotypes, int a, int b);
  double CALPHA(Gene &gene, arma::vec &Y, const std::string &k);
  // Li and Leal 2008
  double CMC(Gene &gene, arma::vec &Y, const std::string &k, double maf = 0.005);
  double CMC1df(Gene &gene, arma::vec &Y, const std::string &k);
  // Morris and Zeggini 2010
  double RVT1(Gene &gene, arma::vec &Y, arma::mat design, arma::vec &initial_beta, const std::string &k, bool linear);
  double RVT2(Gene &gene,
              arma::vec &Y,
              arma::mat design,
              arma::vec &initial_beta,
              const std::string &k,
              bool linear);
  // Wu, Guan, and Pankow 2016
  double SKAT(Gene &gene,
              const std::string &k,
              arma::vec &phenotypes,
              int a,
              int b,
              bool detail,
              bool linear,
              bool permute,
              bool shuffle);
  // Wu, Guan, and Pankow 2016
  double SKATO(Gene &gene,
               const std::string &k,
               arma::vec &phenotypes,
               int a,
               int b,
               bool detail = false,
               bool linear = false);
  double VAAST(Gene &gene,
               arma::vec &Y,
               const std::string &k,
               bool score_only_minor = true,
               bool score_only_alternative = true,
               double site_penalty = 2.0,
               arma::uword group_threshold = 4,
               bool detail = false,
               bool biallelic = false);
  double VT(Gene &gene, const std::string &k, arma::vec &phenotypes);

  // Madsen, Browning 2009
  double WSS(Gene &gene, arma::vec &Y, const std::string &k);

private:
  const std::string method_;

  // SKAT support fields
  Kernel kernel_;
  // Weights for SKAT-O
  static constexpr std::array<double, 8> rho_{0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 1};

  // Check weighting
  void check_weights(Gene &gene, const std::string &k, int a = 1, int b = 25, bool no_weight = false);

  // SKAT Null Model
  std::shared_ptr<SKATR_Null> obj_;
  std::shared_ptr<SKATR_Linear_Null> lin_obj_;

  // VT Helper
  std::shared_ptr<VT_Res> vt_obj_;

  // TaskParams
  TaskParams tp_;
};

#endif //PERMUTE_ASSOCIATE_METHODS_HPP
