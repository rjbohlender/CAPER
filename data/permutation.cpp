//
// Created by Bohlender,Ryan James on 8/4/18.
//

#include <stocc/stocc.h>
#include <ctime>
#include <thread>
#include <cassert>
#include "permutation.hpp"

Permute::Permute()
	: sto(std::random_device{}()) {}

Permute::Permute(int seed)
	: sto(seed) {}

void Permute::get_permutations(std::shared_ptr<std::vector<std::vector<int8_t>>> permutations,
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

std::vector<std::vector<int32_t>> Permute::permutations_maj_bin(int nperm,
																arma::vec &odds,
																arma::uword ncases,
																arma::uvec &mac_indices,
																arma::uvec &maj_indices,
																const std::string &transcript,
																arma::uword maj_nbins,
																double lower_bin_cutoff,
																double upper_bin_cutoff) {
  if (bins_built.find(transcript) == bins_built.end()) {
	bins_built[transcript] = false;
  }

  if (!bins_built[transcript]) {
	ret[transcript] = std::vector<std::vector<int32_t>>(nperm);
	for (int i = 0; i < nperm; i++) {
	  ret[transcript][i] = std::vector<int32_t>(odds.n_rows, 0);
	}
	odds_[transcript] = arma::conv_to<std::vector<double>>::from(odds(mac_indices));
	m[transcript] = std::vector<int32_t>(odds_[transcript].size(), 1);
	sort_mac_idx[transcript] = mac_indices;
	mac_bins[transcript] = mac_indices.n_elem;

	build_major_bins(odds, maj_indices, transcript, maj_nbins, lower_bin_cutoff, upper_bin_cutoff);
  }

  // Generate permutations
  for (int i = 0; i < nperm; i++) {
	std::vector<int32_t> tmp(odds_[transcript].size(), 0);
	sto.MultiFishersNCHyp(&tmp[0], &(m[transcript][0]), &(odds_[transcript][0]), ncases, odds_[transcript].size());

	// Unpack bins
	arma::uword filled = 0;
	for (int j = 0; j < mac_bins[transcript]; j++) { // for each bin
	  arma::uvec
		  r = unpack(tmp[j], m[transcript][j], j < mac_bins[transcript]); // unpack and randomize cases into a vector
	  for (int k = 0; k < r.n_elem; k++) { //
		ret[transcript][i][sort_mac_idx[transcript][k + filled]] = r(k);
	  }
	  filled += m[transcript][j];
	}
	filled = 0;
	for (int j = mac_bins[transcript]; j < mac_bins[transcript] + maj_bins[transcript]; j++) {
	  arma::uvec
		  r = unpack(tmp[j], m[transcript][j], j < mac_bins[transcript]); // Don't have to shuffle maj allele carriers
	  for (int k = 0; k < r.n_elem; k++) {
		ret[transcript][i][sort_maj_idx[transcript][k + filled]] = r(k);
	  }
	  filled += m[transcript][j];
	}
  }

  return ret[transcript];
}

Permute &Permute::operator=(const Permute &rhs) {
  sto = rhs.sto;

  return *this;
}

Permute::Permute(const Permute &other)
	: sto(other.sto) {}

void Permute::permute_thread(std::shared_ptr<std::vector<std::vector<int8_t>>> p,
							 int32_t *m,
							 double *odds,
							 int ncases,
							 int ngroups,
							 int offset,
							 int nperm,
							 int seed) {
  StochasticLib3 rng(seed);

  std::vector<int32_t> tmp(ngroups, 0);

  for (int i = 0; i < nperm; i++) {
	rng.MultiFishersNCHyp(&tmp[0], m, odds, ncases, ngroups);
	for (int j = 0; j < ngroups; j++) {
	  (*p)[offset + i][j] = static_cast<int8_t>(tmp[j]);
	}
  }

}

/**
 * @brief Generate permutations while binning both minor and major allele carriers
 * @param nperm Number of permutations to generate
 * @param odds Odds for all samples
 * @param ncases Total number of cases among samples
 * @param mac_indices Indices of minor allele carriers
 * @param maj_indices Indices of major allele carriers
 * @param transcript Current transcript
 * @param approximate Number of minor allele carrier bins
 * @param maj_nbins Number of major allele carrier bins
 * @param lower_bin_cutoff Independently bin samples with odds below this cutoff
 * @param upper_bin_cutoff Independently bin samples with odds above this cutoff
 * @return
 */
std::vector<std::vector<int32_t>> Permute::permutations_mac_bin(int nperm,
																arma::vec &odds,
																arma::uword ncases,
																arma::uvec &mac_indices,
																arma::uvec &maj_indices,
																const std::string &transcript,
																arma::uword &approximate,
																arma::uword maj_nbins,
																double lower_bin_cutoff,
																double upper_bin_cutoff) {
  if (bins_built.find(transcript) == bins_built.end()) {
	bins_built[transcript] = false;
  }

  if (!bins_built[transcript]) {
	ret[transcript] = std::vector<std::vector<int32_t>>(nperm);
	for (int i = 0; i < nperm; i++) {
	  ret[transcript][i] = std::vector<int32_t>(odds.n_rows, 0);
	}

	// Create minor bins
	arma::uvec odds_sort = arma::sort_index(odds(mac_indices));
	build_minor_bins(odds, mac_indices, transcript, approximate, odds_sort, lower_bin_cutoff, upper_bin_cutoff);

	int msize = std::accumulate(m[transcript].begin(), m[transcript].end(), 0);
	assert(msize == mac_indices.n_elem);

	// Maj bin
	if (maj_indices.n_elem > 0) {
	  build_major_bins(odds, maj_indices, transcript, maj_nbins, lower_bin_cutoff, upper_bin_cutoff);
	}

	msize = std::accumulate(m[transcript].begin(), m[transcript].end(), 0);
	assert(msize == mac_indices.n_elem + maj_indices.n_elem);
	bins_built[transcript] = true;
  }

  // Generate permutations
  for (int i = 0; i < nperm; i++) {
	std::vector<int32_t> tmp(odds_[transcript].size(), 0);
	sto.MultiFishersNCHyp(&tmp[0], &(m[transcript][0]), &(odds_[transcript][0]), ncases, odds_[transcript].size());

	// Unpack bins
	arma::uword filled = 0;
	for (int j = 0; j < mac_bins[transcript]; j++) { // for each bin
	  arma::uvec
		  r = unpack(tmp[j], m[transcript][j], j < mac_bins[transcript]); // unpack and randomize cases into a vector
	  for (int k = 0; k < r.n_elem; k++) { //
		ret[transcript][i][sort_mac_idx[transcript][k + filled]] = r(k);
	  }
	  filled += m[transcript][j];
	}
	filled = 0;
	for (int j = mac_bins[transcript]; j < mac_bins[transcript] + maj_bins[transcript]; j++) {
	  arma::uvec
		  r = unpack(tmp[j], m[transcript][j], j < mac_bins[transcript]); // Don't have to shuffle major allele carriers
	  for (int k = 0; k < r.n_elem; k++) {
		ret[transcript][i][sort_maj_idx[transcript][k + filled]] = r(k);
	  }
	  filled += m[transcript][j];
	}
  }

  return ret[transcript];
}

/**
 * @brief One time construction of bins for major allele carriers for a transcript
 * @param odds The computed odds for all samples
 * @param mac_indices The indices of the minor allele carriers for the current transcript
 * @param transcript The current transcript
 * @param approximate Number of bins to group the minor allele carriers into
 * @param odds_sort The indices of samples in odds sorted order
 * @param lower_bin_cutoff Separate binning for samples with odds lower than the cutoff
 * @param upper_bin_cutoff Separate binning for sample with odds higher than the cutoff
 */
void Permute::build_minor_bins(const arma::vec &odds,
							   const arma::uvec &mac_indices,
							   const std::string &transcript,
							   const arma::uword &approximate,
							   const arma::uvec &odds_sort,
							   double lower_bin_cutoff,
							   double upper_bin_cutoff) {
  sort_mac_idx[transcript] = mac_indices(odds_sort);
  // subset bins
  arma::vec mac_odds = arma::sort(odds(mac_indices));
  arma::vec lower(mac_odds(arma::find(mac_odds < lower_bin_cutoff)));
  arma::vec upper(mac_odds(arma::find(mac_odds > upper_bin_cutoff)));
  arma::vec mac_odds_binnable(mac_odds(arma::find(mac_odds >= lower_bin_cutoff && mac_odds <= upper_bin_cutoff)));

  // Fill in lower
  for (const auto &v : lower) {
	m[transcript].push_back(1);
	odds_[transcript].push_back(v);
  }

  if (mac_odds_binnable.n_elem > 0) {
	double nbins = approximate;
	double bin_width = ((arma::max(mac_odds_binnable) + 0.5) - arma::min(mac_odds_binnable)) / nbins;
	mac_bins[transcript] = nbins + lower.n_elem + upper.n_elem;
	double min_mac_odds = arma::min(mac_odds_binnable);

	for (arma::uword i = 0; i < nbins; i++) {
	  std::vector<arma::uword> odds_in_range;
	  for (arma::uword j = 0; j < mac_odds_binnable.n_elem; j++) { // sorted
		if (mac_odds_binnable(j) >= min_mac_odds + i * bin_width
			&& mac_odds_binnable(j) < min_mac_odds + (i + 1) * bin_width) {
		  odds_in_range.push_back(j);
		} else if (mac_odds_binnable(j) < min_mac_odds + i * bin_width) {
		  continue;
		} else {
		  break;
		}
	  }
	  arma::uvec odds_spanned = arma::conv_to<arma::uvec>::from(odds_in_range);
	  if (odds_spanned.n_elem > 0) {
		// Set number in group
		m[transcript].push_back(odds_spanned.n_elem);
		// Set odds for group
		odds_[transcript].push_back(arma::mean(mac_odds_binnable(odds_spanned)));
	  } else {
		mac_bins[transcript]--;
	  }
	}
  }

  // Fill in upper
  for (const auto &v : upper) {
	m[transcript].push_back(1);
	odds_[transcript].push_back(v);
  }
}

/**
 * @brief One time construction of bins for major allele carriers for a transcript
 * @param odds The computed odds for all samples
 * @param maj_indices The indices of the major allele carriers for the current transcript
 * @param transcript The current transcript
 * @param maj_nbins Number of bins to group the major allele carriers into
 * @param lower_bin_cutoff Separate binning for samples with odds lower than the cutoff
 * @param upper_bin_cutoff Separate binning for sample with odds higher than the cutoff
 */
void Permute::build_major_bins(const arma::vec &odds,
							   const arma::uvec &maj_indices,
							   const std::string &transcript,
							   arma::uword maj_nbins,
							   double lower_bin_cutoff,
							   double upper_bin_cutoff) {
  arma::vec maj_odds = arma::sort(odds(maj_indices));
  arma::vec lower(maj_odds(arma::find(maj_odds < lower_bin_cutoff)));
  arma::vec upper(maj_odds(arma::find(maj_odds > upper_bin_cutoff)));
  arma::vec maj_odds_binnable(maj_odds(arma::find(maj_odds >= lower_bin_cutoff && maj_odds <= upper_bin_cutoff)));

  // Fill in lower
  for (const auto &v : lower) {
	m[transcript].push_back(1);
	odds_[transcript].push_back(v);
  }
  if (maj_odds_binnable.n_elem > 0) {
	double nbins = maj_nbins + lower.n_elem + upper.n_elem;
	maj_bins[transcript] = nbins;
	arma::uword stride = maj_odds_binnable.n_elem / nbins;
	arma::uvec odds_sort = arma::sort_index(maj_odds_binnable);
	sort_maj_idx[transcript] = maj_indices(arma::sort_index(maj_odds));

	for (arma::uword i = 0; i < maj_nbins; i++) {
	  arma::span cur;
	  if (i == maj_nbins - 1) { // Bin remaining, bin may be larger than other bins
		cur = arma::span(i * stride, std::max(i * stride + stride - 1, maj_odds_binnable.n_elem - 1));
	  } else {
		cur = arma::span(i * stride, std::min(i * stride + stride - 1, maj_odds_binnable.n_elem - 1));
	  }
	  arma::vec odds_spanned = maj_odds_binnable(odds_sort(cur));
	  // Set number in group
	  m[transcript].push_back(odds_spanned.n_elem);
	  // Set odds for group
	  odds_[transcript].push_back(arma::mean(odds_spanned));
	}
  }
  // Fill in upper
  for (const auto &v : upper) {
	m[transcript].push_back(1);
	odds_[transcript].push_back(v);
  }
}

/**
 * @brief Unpack the permuted values into random order for a bin
 * @param successes Number of cases within bin
 * @param bin_size Total number of samples in bin
 * @param shuffle Whether to randomize or not
 * @return Vector of phenotype states
 */
auto Permute::unpack(int successes, int bin_size, bool shuffle) -> arma::uvec {
  arma::uvec r(bin_size, arma::fill::zeros);
  if (successes > 0) {
	r(arma::span(0, successes - 1)).fill(1);
  }
  // Fisher-Yates Shuffle
  for (arma::sword i = r.n_elem - 1; i >= 1; --i) {
	auto j = static_cast<arma::sword>(sto.IRandom(0, i));
	arma::uword tmp = r[i];
	r[i] = r[j];
	r[j] = tmp;
  }
  return r;
}

/**
 * @brief Clear the state of the object to free memory
 */
auto Permute::reset() -> void {
  bins_built.clear();
  odds_.clear();
  m.clear();
  mac_bins.clear();
  maj_bins.clear();
  ret.clear();
  sort_mac_idx.clear();
  sort_maj_idx.clear();
}

std::vector<std::vector<int32_t>> Permute::epsilon_permutation(int nperm,
															   arma::vec &odds,
															   arma::uword ncases,
															   const std::string &transcript,
															   double epsilon) {
  if (bins_built.find(transcript) == bins_built.end()) {
	bins_built[transcript] = false;
  }

  if (!bins_built[transcript]) {
	ret[transcript] = std::vector<std::vector<int32_t>>(nperm);
	for (int i = 0; i < nperm; i++) {
	  ret[transcript][i] = std::vector<int32_t>(odds.n_rows, 0);
	}
	sort_idx[transcript] = arma::sort_index(odds);
	arma::vec odds_sorted = odds(sort_idx[transcript]);

	for (double cur = arma::min(odds_sorted);;) {
	  arma::uvec in_range = arma::find(odds_sorted >= cur && odds_sorted < cur + epsilon);
	  m[transcript].push_back(in_range.n_elem);
	  odds_[transcript].push_back(arma::mean(odds_sorted(in_range)));
	  arma::uword next = in_range.max() + 1;
	  if (next < odds_sorted.n_elem) {
		cur = odds_sorted(next);
	  } else {
		break;
	  }
	}
	std::cerr << "nbins: " << m[transcript].size() << std::endl;
	bins_built[transcript] = true;
  }
#ifndef NDEBUG
  arma::uword msize = std::accumulate(m[transcript].begin(), m[transcript].end(), 0);
  assert(msize == odds.n_rows);
#endif
  for (int i = 0; i < nperm; i++) {
	std::vector<int32_t> tmp(odds_[transcript].size(), 0);
	sto.MultiFishersNCHyp(&tmp[0], &(m[transcript][0]), &(odds_[transcript][0]), ncases, odds_[transcript].size());

	// Unpack bins
	arma::uword filled = 0;
	for (int j = 0; j < m[transcript].size(); j++) { // for each bin
	  arma::uvec
		  r = unpack(tmp[j], m[transcript][j], true); // unpack and randomize cases into a vector
	  for (int k = 0; k < r.n_elem; k++) { //
		ret[transcript][i][sort_idx[transcript][k + filled]] = r(k);
	  }
	  filled += m[transcript][j];
	}
  }
  return ret[transcript];
}

