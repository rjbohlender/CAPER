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
  Permute(const Permute &other);
  Permute &operator=(const Permute &rhs);

  void get_permutations(std::shared_ptr<std::vector<std::vector<int32_t>>> permutations,
						arma::colvec &odds,
						arma::uword ncases,
						arma::uword nperm,
						arma::uword nthreads);

  void permute_thread(std::shared_ptr<std::vector<std::vector<int32_t>>> p,
						int32_t *m,
						double *odds,
						int ncases,
						int ngroups,
						int offset,
						int nperm,
						int seed);
  std::vector<std::vector<int32_t>> permutations_maj_bin(int nperm,
														 arma::vec &odds,
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

  StochasticLib3 sto;
};
#endif //PERMUTE_ASSOCIATE_PERMUTATION_HPP
