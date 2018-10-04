//
// Created by Bohlender,Ryan James on 8/4/18.
//

#include <stocc/stocc.h>
#include <ctime>
#include <thread>
#include "permutation.hpp"

Permute::Permute()
	: sto((int32_t) time(nullptr)) {}

void Permute::get_permutations(std::vector<std::vector<int32_t>> *permutations,
							   arma::colvec &odds,
							   int ncases,
							   unsigned long nperm,
							   unsigned long nthreads) {
  std::vector<double> odds_ = arma::conv_to<std::vector<double>>::from(odds);
  std::vector<int32_t> m(odds.n_rows, 1);

  // Initialize permutations
  (*permutations).resize(nperm);
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
									&m[0],
									&odds_[0],
									ncases,
									static_cast<int>(odds.n_rows),
									offset,
									remaining,
									seed));
	} else {
	  threads.push_back(std::thread(&Permute::permute_thread,
									this,
									permutations,
									&m[0],
									&odds_[0],
									ncases,
									static_cast<int>(odds.n_rows),
									offset,
									step,
									seed));
	}
	remaining -= step;
  }

  for (auto &t : threads) {
	t.join();
  }

#if 0
  // Generate permutations
  for (int i = 0; i < nperm; i++) {
	sto.MultiFishersNCHyp(&(*permutations)[i][0], &m[0], &odds_[0], ncases, static_cast<int>(odds.n_rows));
  }
#endif
}

std::vector<std::vector<int32_t>> Permute::cases_in_bins(int nperm,
														 arma::colvec &odds,
														 int ncases,
														 std::vector<int32_t> &bin_counts) {
  int colors = static_cast<int>(odds.n_rows);

  std::vector<double> odds_ = arma::conv_to<std::vector<double>>::from(odds);
  std::vector<int32_t> m(odds.n_rows);
  std::vector<std::vector<int32_t>> ret(nperm);

  for (int i = 0; i < odds.n_rows; i++) {
	m[i] = bin_counts[i];
  }

  for (int i = 0; i < nperm; i++) {
	ret[i] = std::vector<int32_t>(odds.n_rows);
  }

  // Generate permutations
  for (int i = 0; i < nperm; i++) {
	sto.MultiFishersNCHyp(&ret[i][0], &m[0], &odds_[0], ncases, colors);
  }

  return ret;
}

std::vector<int32_t> Permute::random_case_count(int nperm,
												arma::uvec &mac_indices,
												arma::uvec &maj_indices,
												arma::vec &prob,
												int ncases) {
  std::vector<int32_t> ret = std::vector<int32_t>(nperm);

  // bin1 for minor allele carriers
  // bin2 for major allele carriers
  int32_t bin1_count = static_cast<int32_t>(mac_indices.size());
  int32_t bin2_count = static_cast<int32_t>(maj_indices.size());

  double bin1_odds = arma::mean(prob(mac_indices)) / (1 - arma::mean(prob(mac_indices)));
  double bin2_odds = arma::mean(prob(maj_indices)) / (1 - arma::mean(prob(maj_indices)));

  int32_t sample[2] = {0, 0};
  int32_t m[2] = {bin1_count, bin2_count};
  double odds_[2] = {bin1_odds, bin2_odds};

  for (int i = 0; i < nperm; i++) {
	sto.MultiFishersNCHyp(sample, m, odds_, ncases, 2);
	// We only return the case counts among minor allele carriers
	ret[i] = sample[0];
  }
  return ret;
}

std::vector<int32_t> Permute::random_case_count(int nperm,
												arma::uvec &mac_indices,
												arma::uvec &maj_indices,
												arma::vec &prob,
												int ncases,
												int n_maj_bins) {
  std::vector<int32_t> sample(n_maj_bins + 1, 0);
  std::vector<int32_t> m(n_maj_bins + 1);
  std::vector<double> odds(n_maj_bins + 1, 0);
  std::vector<int32_t> ret = std::vector<int32_t>(nperm);
  std::vector<arma::span> span_vec;
  arma::vec maj_prob = prob(maj_indices);
  arma::uword split = maj_indices.n_rows / n_maj_bins;

  // Set splits
  for (int i = 0; i < n_maj_bins; i++) {
	if (i == n_maj_bins - 1) {
	  span_vec.emplace_back(arma::span(i * split, maj_indices.n_rows - 1));
	} else {
	  span_vec.emplace_back(arma::span(i * split, (i + 1) * split - 1));
	}
  }

  // Set bin counts and odds
  m[0] = mac_indices.size();
  odds[0] = arma::mean(prob(mac_indices)) / (1 - arma::mean(prob(mac_indices)));
  for (int i = 0; i < m.size() - 1; i++) {
	m[i + 1] = span_vec[i].b - span_vec[i].a;
	odds[i + 1] = arma::mean(maj_prob(span_vec[i])) / (1 - arma::mean(maj_prob(span_vec[i])));
  }

  for (int i = 0; i < nperm; i++) {
	sto.MultiFishersNCHyp(&sample[0], &m[0], &odds[0], ncases, sample.size());
	// We only return the case counts among minor allele carriers
	ret[i] = sample[0];
  }
  return ret;
}

arma::vec Permute::calculate_fisher_mean(int32_t n, arma::vec &odds) {
  int colors = static_cast<int>(odds.n_rows); // One color each
  std::vector<int32_t> m(colors, 1);
  std::vector<double> mu(colors, 0);
  std::vector<double> odds_ = arma::conv_to<std::vector<double>>::from(odds);

  CMultiFishersNCHypergeometric mfnch(n, &m[0], &odds_[0], colors);

  mfnch.mean(&mu[0]);

  return arma::conv_to<arma::vec>::from(mu);
}

std::vector<std::vector<int32_t>> Permute::permutations_maj_bin(int nperm,
																arma::vec &prob,
																int ncases,
																arma::uvec &mac_indices,
																arma::uvec &maj_indices) {
  std::vector<std::vector<int32_t>> ret(nperm);
  for (int i = 0; i < nperm; i++) {
	ret[i] = std::vector<int32_t>(prob.n_rows);
  }

  std::vector<double> odds = arma::conv_to<std::vector<double>>::from(prob(mac_indices) / (1 - prob(mac_indices)));
  arma::vec maj_prob = arma::sort(prob(maj_indices));

  // Create bins
  double nbins = 1.;
  std::vector<int32_t> m(odds.size(), 1);
  for (double i = 0; i < nbins; i++) {
	arma::uvec cur = arma::find(maj_prob >= i / nbins && maj_prob < (i + 1.) / nbins);
	// Set number in group
	m.push_back(cur.size());
	// Set odds for group
	if (cur.size() == 0) {
	  odds.push_back(0);
	} else {
	  odds.push_back(arma::mean(maj_prob(cur) / (1 - maj_prob(cur))));
	}
  }

  // Generate permutations
  for (int i = 0; i < nperm; i++) {
	sto.MultiFishersNCHyp(&ret[i][0], &m[0], &odds[0], ncases, odds.size());
  }

  return ret;
}

Permute &Permute::operator=(const Permute &rhs) {
  sto = rhs.sto;

  return *this;
}

Permute::Permute(const Permute &other)
	: sto(other.sto) {}

void Permute::permute_thread(std::vector<std::vector<int32_t>> *p,
							 int32_t *m,
							 double *odds,
							 int ncases,
							 int ngroups,
							 int offset,
							 int nperm,
							 int seed) {
  StochasticLib3 rng(seed);

  for (int i = 0; i < nperm; i++) {
	rng.MultiFishersNCHyp(&(*p)[offset + i][0], m, odds, ncases, ngroups);
  }

}

