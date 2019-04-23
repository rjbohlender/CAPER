//
// Created by Bohlender,Ryan James on 10/2/18.
//

#ifndef PERMUTE_ASSOCIATE_VAAST_HPP
#define PERMUTE_ASSOCIATE_VAAST_HPP

#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

#include "../data/gene.hpp"
#include "../data/covariates.hpp"

class VariantGroup {
public:
  arma::uword nvariants;
  arma::sp_mat Xcollapse;
  arma::vec Y;
  arma::vec weight;
  const double site_penalty;
  const bool som;
  const bool soa;

  arma::vec case_allele1;
  arma::vec control_allele1;
  arma::vec case_allele0;
  arma::vec control_allele0;

  arma::uword n_case;
  arma::uword n_control;

  arma::vec variant_scores;
  arma::uvec mask;

  double score;

  VariantGroup(arma::sp_mat X,
				 arma::vec Y,
				 arma::vec weights,
				 arma::uword group_threshold,
				 double site_penalty,
				 bool score_only_minor,
				 bool score_only_alternative);

private:
  void variant_mask();

  double Score(const arma::sp_mat &X, const arma::vec &Y, const arma::vec &w);
  arma::vec LRT();
  arma::vec log_likelihood(arma::vec &freq, arma::vec &allele0, arma::vec &allele1);
};

struct VAAST {
  const bool som;      // score_only_minor
  const bool soa;      // score_only_alternative
  const bool detail;   // Add detailed output
  const bool biallelic; // Biallelic variants get additional score
  const std::string k; // Transcript
  const arma::uword group_threshold;
  const double site_penalty;

  arma::sp_mat X;
  const arma::vec Y;

  arma::vec weights;

  arma::vec case_allele1;
  arma::vec control_allele1;
  arma::vec case_allele0;
  arma::vec control_allele0;

  arma::uword n_case;
  arma::uword n_control;

  arma::uvec mask;

  std::vector<VariantGroup> groups;
  arma::vec vaast_site_scores;
  arma::vec expanded_scores;

  double score;

  // Constructors
  VAAST(Gene &gene,
		arma::vec &Y,
		const std::string &k,
		bool score_only_minor,
		bool score_only_alternative,
		double site_penalty,
		arma::uword group_threshold,
		bool detail,
		bool biallelic);
  VAAST(arma::sp_mat X,
		arma::vec &Y,
		arma::vec &weights,
		std::vector<std::string> &positions_,
		const std::string &k,
		bool score_only_minor,
		bool score_only_alternative,
		bool biallelic,
		arma::uword group_threshold,
		double site_penalty);

  void check_weights(Gene &gene);
  double Score();
  double Score(const arma::sp_mat &X, const arma::vec &Y, const arma::vec &w);
  arma::vec LRT();
  arma::vec log_likelihood(arma::vec &freq, arma::vec &allele0, arma::vec &allele1);
  void variant_grouping(const arma::sp_mat &X,
						  const arma::vec &Y,
						  const arma::vec &w,
						  std::vector<std::string> &positions);
  void variant_bitmask(const arma::sp_mat &X, const arma::vec &Y, const arma::vec &w);

  arma::uvec setdiff(arma::uvec x, arma::uvec y);
};

#endif //PERMUTE_ASSOCIATE_VAAST_HPP
