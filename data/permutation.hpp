//
// Created by Bohlender,Ryan James on 8/4/18.
//

#ifndef PERMUTE_ASSOCIATE_PERMUTATION_HPP
#define PERMUTE_ASSOCIATE_PERMUTATION_HPP

#include <armadillo>
#include <stocc/randomc.h>
#include <stocc/stocc.h>

struct Permute {
  Permute();

  std::vector<std::vector<int32_t>> get_permutations(unsigned long nperm, arma::colvec &odds, int ncases);
  std::vector<std::vector<int32_t>> permutations_maj_bin(int nperm,
														 arma::vec &prob,
														 int ncases,
														 arma::uvec &mac_indices,
														 arma::uvec &maj_indices);
  arma::vec calculate_fisher_mean(int32_t n, arma::vec &odds);
  std::vector<std::vector<int32_t>> cases_in_bins(int nperm,
												  arma::colvec &odds,
												  int ncases,
												  std::vector<int32_t> &bin_counts);
  std::vector<int32_t> random_case_count(int nperm,
										 arma::uvec &mac_indices,
										 arma::uvec &maj_indices,
										 arma::vec &prob,
										 int ncases);
  std::vector<int32_t> random_case_count(int nperm,
										 arma::uvec &mac_indices,
										 arma::uvec &maj_indices,
										 arma::vec &prob,
										 int ncases,
										 int n_maj_bins);

  StochasticLib3 sto_rcc;
  StochasticLib3 sto_cib;
  StochasticLib3 sto_perm;
  StochasticLib3 sto_pmj;
};
#endif //PERMUTE_ASSOCIATE_PERMUTATION_HPP
