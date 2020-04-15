//
// Created by Bohlender,Ryan James on 8/4/18.
//

#ifndef PERMUTE_ASSOCIATE_PERMUTATION_HPP
#define PERMUTE_ASSOCIATE_PERMUTATION_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <stocc/randomc.h>
#include <stocc/stocc.h>
#include <memory>


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
														 arma::uword ncases,
														 arma::uvec &mac_indices,
														 arma::uvec &maj_indices,
														 const std::string &transcript,
														 arma::uword maj_nbins,
														 double lower_bin_cutoff,
														 double upper_bin_cutoff);
  std::vector<std::vector<int32_t>> permutations_mac_bin(int nperm,
														 arma::vec &odds,
														 arma::uword ncases,
														 arma::uvec &mac_indices,
														 arma::uvec &maj_indices,
														 const std::string &transcript,
														 arma::uword &approximate,
														 arma::uword maj_nbins,
														 double lower_bin_cutoff,
														 double upper_bin_cutoff);

  auto unpack(int successes, int bin_size, bool shuffle) -> arma::uvec;
  auto reset() -> void;

  StochasticLib3 sto;
  // Preserve group info for transcript
  std::map<std::string, bool> bins_built; // Whether the bins have been built or not -- for each transcript
  std::map<std::string, std::vector<double>> odds_; // Odds for the bins -- averaged within bin
  std::map<std::string, std::vector<int32_t>> m; // Number of samples in each bin
  std::map<std::string, double> mac_bins; // Total number of minor allele carrier bins = approximate + lower + upper or number of minor allele carriers
  std::map<std::string, double> maj_bins; // Total number of major allele carrier bins = maj_nbins + lower + upper
  std::map<std::string, std::vector<std::vector<int32_t>>> ret;
  std::map<std::string, arma::uvec> sort_mac_idx; // Indices of minor allele carriers in odds sorted order
  std::map<std::string, arma::uvec> sort_maj_idx; // Indices of major allele carriers in odds sorted order
  std::map<std::string, std::vector<arma::uvec>> mac_spans;
  std::map<std::string, std::vector<arma::span>> maj_spans;
  void build_major_bins(const arma::vec &odds,
						const arma::uvec &maj_indices,
						const std::string &transcript,
						arma::uword maj_nbins,
						double lower_bin_cutoff,
						double upper_bin_cutoff);
  void build_minor_bins(const arma::vec &odds,
						const arma::uvec &mac_indices,
						const std::string &transcript,
						const arma::uword &approximate,
						const arma::uvec &odds_sort,
						double lower_bin_cutoff,
						double upper_bin_cutoff);
};
#endif //PERMUTE_ASSOCIATE_PERMUTATION_HPP
