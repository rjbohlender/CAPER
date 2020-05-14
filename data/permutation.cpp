//
// Created by Bohlender,Ryan James on 8/4/18.
//

#include <stocc/stocc.h>
#include <ctime>
#include <cassert>
#include "permutation.hpp"

Permute::Permute()
	: sto(std::random_device{}()), bins_built(false)  {}

Permute::Permute(int seed)
	: sto(seed), bins_built(false) {}

void Permute::generate_permutations(std::shared_ptr<std::vector<std::vector<int8_t>>> permutations,
									arma::colvec &odds,
									arma::uword ncases,
									arma::uword nperm,
									arma::uword nthreads,
									double epsilon) {
  if (!bins_built) {
	sort_idx = arma::sort_index(odds);
	arma::vec odds_sorted = odds(sort_idx);

	for (double cur = arma::min(odds_sorted);;) {
	  arma::uvec
		  in_range = arma::find(odds_sorted >= cur && odds_sorted < cur + epsilon + std::numeric_limits<double>::min());
	  m.push_back(in_range.n_elem);
	  bin_odds.push_back(arma::mean(odds_sorted(in_range)));
	  arma::uword next = in_range.max() + 1;
	  if (next < odds_sorted.n_elem) {
		cur = odds_sorted(next);
	  } else {
		break;
	  }
	}
	std::cerr << "nbins: " << m.size() << std::endl;
	nsamples = sort_idx.n_elem;
	bins_built = true;
  }
#ifndef NDEBUG
  arma::uword msize = std::accumulate(m.begin(), m.end(), 0);
  assert(msize == odds.n_rows);
#endif

  // Initialize permutations
  permutations->resize(nperm);
  for (int i = 0; i < nperm; i++) {
	(*permutations)[i].resize(odds.n_rows);
  }

  int step = nperm / nthreads;
  int remaining = nperm;
  std::vector<std::thread> threads;
  for (int i = 0; i < nthreads; i++) {
	int seed = sto.IRandom(0, std::numeric_limits<int>::max());
	int offset = i * step;
	if (remaining < 0) {
	  std::cerr << "Failed during permutation.\n";
	  std::exit(-1);
	}
	if (i == nthreads - 1) {
	  threads.push_back(std::thread(&Permute::permute_thread,
									this,
									permutations,
									ncases,
									offset,
									remaining,
									seed));
	} else {
	  threads.push_back(std::thread(&Permute::permute_thread,
									this,
									permutations,
									ncases,
									offset,
									step,
									seed));
	}
	remaining -= step;
  }

  for (auto &t : threads) {
	t.join();
  }
}

void Permute::permute_thread(std::shared_ptr<std::vector<std::vector<int8_t>>> p,
							 int ncases,
							 int offset,
							 int nperm,
							 int seed) {
  StochasticLib3 rng(seed);
  std::vector<int32_t> tmp(bin_odds.size(), 0);

  for (int i = 0; i < nperm; i++) {
	rng.MultiFishersNCHyp(&tmp[0], &(m[0]), &(bin_odds[0]), ncases, bin_odds.size());
	// Unpack bins
	arma::uword filled = 0;
	for (int j = 0; j < m.size(); j++) { // for each bin
	  std::vector<int8_t> r = unpack(tmp[j], m[j], true, rng); // unpack and randomize cases into a vector
	  for (int k = 0; k < r.size(); k++) {
		(*p)[offset + i][sort_idx[k + filled]] = r[k];
	  }
	  filled += m[j];
	}
  }
}

/**
 * @brief Unpack the permuted values into random order for a bin
 * @param successes Number of cases within bin
 * @param bin_size Total number of samples in bin
 * @param shuffle Whether to randomize or not
 * @return Vector of phenotype states
 */
std::vector<int8_t> Permute::unpack(int successes, int bin_size, bool shuffle, StochasticLib3 &rng) {
  std::vector<int8_t> r(bin_size, 0);
  if (successes > 0) {
	for (int i = 0; i < successes; i++) {
	  r[i] = 1;
	}
  }
  // Fisher-Yates Shuffle
  if (shuffle) {
	for (int i = r.size() - 1; i >= 1; --i) {
	  auto j = rng.IRandom(0, i);
	  arma::uword tmp = r[i];
	  r[i] = r[j];
	  r[j] = tmp;
	}
  }
  return r;
}

std::vector<std::vector<int8_t>> Permute::epsilon_permutation(int nperm,
															  arma::vec &odds,
															  arma::uword ncases,
															  const std::string &transcript,
															  double epsilon) {
  if (!bins_built) {
	ret = std::vector<std::vector<int8_t>>(nperm);
	for (int i = 0; i < nperm; i++) {
	  ret[i] = std::vector<int8_t>(odds.n_rows, 0);
	}
	sort_idx = arma::sort_index(odds);
	arma::vec odds_sorted = odds(sort_idx);

	for (double cur = arma::min(odds_sorted);;) {
	  arma::uvec in_range = arma::find(odds_sorted >= cur && odds_sorted < cur + epsilon);
	  m.push_back(in_range.n_elem);
	  bin_odds.push_back(arma::mean(odds_sorted(in_range)));
	  arma::uword next = in_range.max() + 1;
	  if (next < odds_sorted.n_elem) {
		cur = odds_sorted(next);
	  } else {
		break;
	  }
	}
	std::cerr << "nbins: " << m.size() << std::endl;
	bins_built = true;
  }
#ifndef NDEBUG
  arma::uword msize = std::accumulate(m.begin(), m.end(), 0);
  assert(msize == odds.n_rows);
#endif
  for (int i = 0; i < nperm; i++) {
	std::vector<int32_t> tmp(bin_odds.size(), 0);
	sto.MultiFishersNCHyp(&tmp[0], &(m[0]), &(bin_odds[0]), ncases, bin_odds.size());

	// Unpack bins
	arma::uword filled = 0;
	for (int j = 0; j < m.size(); j++) { // for each bin
	  std::vector<int8_t> r = unpack(tmp[j], m[j], true, sto); // unpack and randomize cases into a vector
	  for (int k = 0; k < r.size(); k++) { //
		ret[i][sort_idx[k + filled]] = r[k];
	  }
	  filled += m[j];
	}
  }
  return ret;
}

void Permute::fisher_yates(std::vector<int8_t> &x, StochasticLib3 &rng) {
  for (int i = x.size() - 1; i >= 1; --i) {
	auto j = rng.IRandom(0, i);
	arma::uword tmp = x[i];
	x[i] = x[j];
	x[j] = tmp;
  }
}

