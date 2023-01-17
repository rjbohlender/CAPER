//
// Created by Bohlender,Ryan James on 10/2/18.
//

#ifndef PERMUTE_ASSOCIATE_VAAST_HPP
#define PERMUTE_ASSOCIATE_VAAST_HPP

#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

#include "../data/gene.hpp"
#include "../data/covariates.hpp"

struct Variant {
  double case_allele1 = 0;
  double case_allele0 = 0;
  double control_allele1 = 0;
  double control_allele0 = 0;
  double weight = 1;
  std::string type = "";
  std::string loc = "";
  int index = 0;
  double soft_maf_filter = 0.5;

  // Calculated here
  double score = 0;

  Variant(std::string type, double soft_maf_filter);
  Variant(double case_allele1,
		  double case_allele0,
		  double control_allele1,
		  double control_allele0,
		  double weight,
		  std::string type,
		  std::string loc,
		  int index,
		  double soft_maf_filter);

  void merge(Variant &other);
  void calc_score();
};

struct VAASTLogic {
  const bool detail;   // Add detailed output
  const bool biallelic; // Biallelic variants get additional score
  const std::string k; // Transcript
  const arma::uword group_threshold;
  const double site_penalty;
  bool printed_mergeinfo = false;
  double soft_maf_filter;
  bool legacy = false;

  arma::sp_mat X;
  const arma::vec Y;

  arma::vec weights;

  arma::vec case_allele1;
  arma::vec control_allele1;
  arma::vec case_allele0;
  arma::vec control_allele0;

  arma::uword n_case;
  arma::uword n_control;

  arma::vec vaast_site_scores;
  arma::vec expanded_scores;

  double score;

  // Constructors
  VAASTLogic(Gene &gene,
			 arma::vec &Y_,
			 const std::string &ts,
			 double site_penalty,
			 arma::uword group_threshold,
			 bool detail,
			 bool biallelic,
			 double smf,
			 bool legacy_);
  VAASTLogic(arma::sp_mat X_,
			 arma::vec &Y_,
			 arma::vec &weights,
			 std::vector<std::string> &positions_,
			 std::string k,
			 bool biallelic,
			 arma::uword group_threshold,
			 double site_penalty,
			 double smf,
			 bool legacy_);

  double Score();
  double Score(const arma::sp_mat &X, const arma::vec &Y, const arma::vec &w);
  arma::vec LRT();
  arma::vec log_likelihood(arma::vec &freq, arma::vec &allele0, arma::vec &allele1);
  void alternate_grouping(const arma::sp_mat &X,
					   const arma::vec &Y,
					   const arma::vec &w,
					   std::vector<std::string> &positions);
  void vaast_grouping(const arma::sp_mat &X,
					   const arma::vec &Y,
					   const arma::vec &w,
					   std::vector<std::string> &positions);
};

#endif //PERMUTE_ASSOCIATE_VAAST_HPP
