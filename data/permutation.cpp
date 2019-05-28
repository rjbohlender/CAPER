//
// Created by Bohlender,Ryan James on 8/4/18.
//

#include <stocc/stocc.h>
#include <ctime>
#include <thread>
#include <cassert>
#include "permutation.hpp"

Permute::Permute()
	: sto((int32_t) time(nullptr)) {}

void Permute::get_permutations(std::shared_ptr<std::vector<std::vector<int32_t>>> permutations,
							   arma::colvec &odds,
							   arma::uword ncases,
							   arma::uword nperm,
							   arma::uword nthreads) {
  std::vector<double> odds_ = arma::conv_to<std::vector<double>>::from(odds);
  std::vector<int32_t> m(odds.n_rows, 1);

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
																arma::vec &odds,
																arma::uword ncases,
																arma::uvec &mac_indices,
																arma::uvec &maj_indices,
																const std::string &transcript) {
  if (bins_built.find(transcript) == bins_built.end()) {
	bins_built[transcript] = false;
  }

  if(!bins_built[transcript]) {
	ret[transcript] = std::vector<std::vector<int32_t>>(nperm);
	for (int i = 0; i < nperm; i++) {
	  ret[transcript][i] = std::vector<int32_t>(odds.n_rows, 0);
	}

	odds_[transcript] = arma::conv_to<std::vector<double>>::from(odds(mac_indices));
	arma::uvec odds_sort = arma::sort_index(odds(maj_indices));

	// Create bins
	double nbins = 1;
	arma::uword stride = maj_indices.n_elem / nbins;
	m[transcript] = std::vector<int32_t>(odds_[transcript].size(), 1);
	for (arma::uword i = 0; i < nbins; i++) {
	  arma::span cur(i * stride, std::min(i * stride + stride - 1, maj_indices.n_elem));
	  arma::vec odds_spanned = odds(odds_sort(cur));
	  // Set number in group
	  m[transcript].push_back(odds_spanned.n_elem);
	  // Set odds for group
	  odds_[transcript].push_back(arma::mean(odds_spanned));
	}
	bins_built[transcript] = true;
  }

  // Generate permutations
  for (int i = 0; i < nperm; i++) {
	sto.MultiFishersNCHyp(&(ret[transcript][i][0]), &(m[transcript][0]), &(odds_[transcript][0]), ncases, odds_[transcript].size());
  }

  return ret[transcript];
}

Permute &Permute::operator=(const Permute &rhs) {
  sto = rhs.sto;

  return *this;
}

Permute::Permute(const Permute &other)
	: sto(other.sto) {}

void Permute::permute_thread(std::shared_ptr<std::vector<std::vector<int32_t>>> p,
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
std::vector<std::vector<int32_t>> Permute::permutations_mac_bin(int nperm,
																arma::vec &odds,
																arma::uword ncases,
																arma::uvec &mac_indices,
																arma::uvec &maj_indices,
																arma::uword &approximate,
																const std::string &transcript) {
  if(bins_built.find(transcript) == bins_built.end()) {
	bins_built[transcript] = false;
  }

  if(!bins_built[transcript]) {
	ret[transcript] = std::vector<std::vector<int32_t>>(nperm);
	for (int i = 0; i < nperm; i++) {
	  ret[transcript][i] = std::vector<int32_t>(odds.n_rows, 0);
	}

	// Create minor bins
	arma::uvec odds_sort = arma::sort_index(odds(mac_indices));
   // subset bins
	arma::vec mac_odds = arma::sort(odds(mac_indices));
	double nbins = approximate;
	double bin_width = (arma::max(mac_odds) - arma::min(mac_odds)) / nbins;
	final_bins[transcript] = nbins;
	double min_mac_odds = arma::min(mac_odds);
	for (arma::uword i = 0; i < nbins; i++) {
	  arma::uvec odds_lower = arma::find(mac_odds >= min_mac_odds + i * bin_width);
	  arma::uvec odds_upper = arma::find(mac_odds <= min_mac_odds + (i + 1) * bin_width);
	  arma::uvec odds_spanned = arma::intersect(odds_lower, odds_upper);
	  if(odds_spanned.n_elem > 0) {
		// Set number in group
		m[transcript].push_back(odds_spanned.n_elem);
		// Set odds for group
		odds_[transcript].push_back(arma::mean(mac_odds(odds_spanned)));
	  } else {
	    final_bins[transcript]--;
	  }
	}
	assert(final_bins[transcript] == odds_[transcript].size() && final_bins[transcript] == m[transcript].size());
	// Maj bin
	nbins = 1;
	arma::uword stride = maj_indices.n_elem / nbins;
	odds_sort = arma::sort_index(odds(maj_indices));
	for (arma::uword i = 0; i < nbins; i++) {
	  arma::span cur(i * stride, std::min(i * stride + stride - 1, maj_indices.n_elem));
	  arma::vec odds_spanned = odds(odds_sort(cur));
	  // Set number in group
	  m[transcript].push_back(odds_spanned.n_elem);
	  // Set odds for group
	  odds_[transcript].push_back(arma::mean(odds_spanned));
	}
	bins_built[transcript] = true;
  }
  // successes, bin_size
  auto unpack = [this](int n, int m){
	arma::uvec r(m, arma::fill::zeros);
	if (n > 0) {
	  r(arma::span(0, n - 1)).fill(1);
	}
	// Fisher-Yates Shuffle
	for(arma::sword i = r.n_elem - 1; i >= 1; --i) {
	  auto j = static_cast<arma::sword>(sto.IRandom(0, i));
	  arma::uword tmp = r[i];
	  r[i] = r[j];
	  r[j] = tmp;
	}
	return r;
  };

  // Generate permutations
  for (int i = 0; i < nperm; i++) {
    std::vector<int32_t> tmp(odds_[transcript].size(), 0);
	sto.MultiFishersNCHyp(&tmp[0], &(m[transcript][0]), &(odds_[transcript][0]), ncases, odds_[transcript].size());

	assert(tmp.size() == odds_[transcript].size() && tmp.size() == m[transcript].size());
	// Unpack bins
	arma::uword filled = 0;
	for (int j = 0; j < final_bins[transcript]; j++) {
	  arma::uvec r = unpack(tmp[j], m[transcript][j]);
	  for (int k = 0; k < r.n_elem; k++) {
		ret[transcript][i][k + filled] = r(k);
	  }
	  filled += m[transcript][j];
	}
	ret[transcript][i][filled] = tmp.back();
  }

  return ret[transcript];
}

std::vector<std::vector<int32_t>> Permute::permutations_bin(int nperm,
															arma::vec &odds,
															arma::uword ncases,
															arma::uvec &mac_indices,
															arma::uvec &maj_indices,
															arma::uword &approximate,
															const std::string &transcript) {
  if(bins_built.find(transcript) == bins_built.end()) {
    bins_built[transcript] = false;
  }
  if(!bins_built[transcript]) {
	ret[transcript] = std::vector<std::vector<int32_t>>(nperm);
	for (int i = 0; i < nperm; i++) {
	  ret[transcript][i] = std::vector<int32_t>(odds.n_rows, 0);
	}

	// Create minor bins
	arma::uvec odds_sort = arma::sort_index(odds);
	// subset bins
	double nbins = approximate;
	double bin_width = (arma::max(odds) - arma::min(odds)) / nbins;
	final_bins[transcript] = nbins;
	double min_odds = arma::min(odds);
	for (arma::uword i = 0; i < nbins; i++) {
	  arma::uvec odds_lower = arma::find(odds >= min_odds + i * bin_width);
	  arma::uvec odds_upper = arma::find(odds <= min_odds + (i + 1) * bin_width);
	  arma::uvec odds_spanned = arma::intersect(odds_lower, odds_upper);
	  if(odds_spanned.n_elem > 0) {
		// Set number in group
		m[transcript].push_back(odds_spanned.n_elem);
		// Set odds for group
		odds_[transcript].push_back(arma::mean(odds(odds_spanned)));
	  } else {
		final_bins[transcript]--;
	  }
	}
	bins_built[transcript] = true;
  }

  // successes, bin_size
  auto unpack = [this](int n, int m){
	arma::uvec r(m, arma::fill::zeros);
	if (n > 0) {
	  r(arma::span(0, n - 1)).fill(1);
	}
	// Fisher-Yates Shuffle
	for(arma::sword i = r.n_elem - 1; i >= 1; i--) {
	  auto j = static_cast<arma::sword>(sto.IRandom(0, i));
	  arma::uword tmp = r[i];
	  r[i] = r[j];
	  r[j] = tmp;
	}
	return r;
  };

  for (int i = 0; i < nperm; i++) {
	std::vector<int32_t> tmp(odds_[transcript].size(), 0);
	sto.MultiFishersNCHyp(&tmp[0], &(m[transcript][0]), &(odds_[transcript][0]), ncases, odds_[transcript].size());

	// Unpack bins
	arma::uword filled = 0;
	for (int j = 0; j < final_bins[transcript]; j++) {
	  arma::uvec r = unpack(tmp[j], m[transcript][j]);
	  for (int k = 0; k < r.n_elem; k++) {
		  ret[transcript][i][k + filled] = r(k);
	  }
	  filled += m[transcript][j];
	}
  }
  return ret[transcript];
}

