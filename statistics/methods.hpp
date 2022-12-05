//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_METHODS_HPP
#define PERMUTE_ASSOCIATE_METHODS_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <map>
#include <memory>
#include <string>
#include <tuple>

#include "../data/covariates.hpp"
#include "../data/gene.hpp"
#include "skat.hpp"
#include "vt.hpp"

arma::vec rank(arma::vec &v, const char *direction);

enum class Kernel { Linear, wLinear };

class Methods {
public:
  Methods(const TaskParams &tp, const std::shared_ptr<Covariates> &cov);
  Methods(const TaskParams &tp, const Covariates &cov);

  double call(Gene &gene, Covariates &cov, arma::vec &phenotypes,
              const std::string &transcript, bool detail);

  std::string str();

  void clear(std::vector<std::string> &v);

  // Wu, Guan, and Pankow 2016
  double BURDEN(Gene &gene, const std::string &ts, arma::vec &phenotypes);
  static double CALPHA(Gene &gene, arma::vec &Y, const std::string &ts);
  // Li and Leal 2008
  double CMC(Gene &gene, arma::vec &Y, const std::string &ts,
             double maf = 0.005) const;
  double CMC1df(Gene &gene, arma::vec &Y, const std::string &ts) const;
  // Morris and Zeggini 2010
  double RVT1(Gene &gene, arma::vec &Y, arma::mat design,
              arma::vec &initial_beta, const std::string &ts, bool linear);
  static double RVT2(Gene &gene, arma::vec &Y, arma::mat design,
                     arma::vec &initial_beta, const std::string &ts,
                     bool linear);
  // Wu, Guan, and Pankow 2016
  double SKAT(Gene &gene, const std::string &transcript, arma::vec &phenotypes,
              int a, int b, bool detail, bool linear, bool permute);
  // Wu, Guan, and Pankow 2016
  double SKATO(Gene &gene, const std::string &transcript, arma::vec &phenotypes,
               int a, int b, bool detail = false, bool linear = false);
  double VAAST(Gene &gene, arma::vec &Y, const std::string &k,
               double site_penalty, arma::uword group_threshold, bool detail,
               bool biallelic, double control_freq_cutoff, bool legacy);
  double VT(Gene &gene, const std::string &ts, arma::vec &phenotypes);

  // Madsen, Browning 2009
  static double WSS(Gene &gene, arma::vec &Y, const std::string &k);

private:
  const std::string method_;

  // SKAT support fields
  Kernel kernel_;
  // Weights for SKAT-O
  static constexpr std::array<double, 8> rho_{0,    0.01, 0.04, 0.09,
                                              0.16, 0.25, 0.5,  1};

  // Check weighting
  void check_weights(Gene &gene, const std::string &transcript, int a = 1,
                     int b = 25, bool no_weight = false);

  // SKAT Null Model
  std::shared_ptr<SKATR_Null> obj_;
  std::shared_ptr<SKATR_Linear_Null> lin_obj_;

  // VT Helper
  std::shared_ptr<VT_Res> vt_obj_;

  // TaskParams
  const TaskParams &tp;
};

#endif // PERMUTE_ASSOCIATE_METHODS_HPP
