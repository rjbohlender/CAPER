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
#include <thread>

struct Permute {
  Permute();
  explicit Permute(int seed);

  void generate_permutations(std::shared_ptr<std::vector<std::vector<int8_t>>> permutations,
							 arma::colvec &odds,
							 arma::uword ncases,
							 arma::uword nperm,
							 arma::uword nthreads,
							 double epsilon);

  void permute_thread(std::shared_ptr<std::vector<std::vector<int8_t>>> p,
					  int ncases,
					  int offset,
					  int nperm,
					  int seed);
  std::vector<std::vector<int8_t>> epsilon_permutation(int nperm,
													   arma::vec &odds,
													   arma::uword ncases,
													   const std::string &transcript,
													   double epsilon = 0.01);

  std::vector<int8_t> unpack(int successes, int bin_size, bool shuffle, StochasticLib3 &rng);
  static void fisher_yates(std::vector<int8_t> &x, StochasticLib3 &rng);

  StochasticLib3 sto;
  // Preserve group info for transcript
  bool bins_built; // Whether the bins have been built or not -- for each transcript
  std::vector<double> bin_odds; // Odds for the bins -- averaged within bin
  std::vector<int32_t> m; // Number of samples in each bin
  std::vector<std::vector<int8_t>> ret;
  arma::uvec sort_idx; // Indices of all samples in odds sorted order
  unsigned long long nsamples;
};
#endif //PERMUTE_ASSOCIATE_PERMUTATION_HPP
