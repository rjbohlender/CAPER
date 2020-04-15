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
																const std::string &transcript,
																arma::uword maj_nbins,
																double lower_bin_cutoff,
																double upper_bin_cutoff) {
  if (bins_built.find(transcript) == bins_built.end()) {
	bins_built[transcript] = false;
  }

  if(!bins_built[transcript]) {
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
	  arma::uvec r = unpack(tmp[j], m[transcript][j], j < mac_bins[transcript]); // unpack and randomize cases into a vector
	  for (int k = 0; k < r.n_elem; k++) { //
		ret[transcript][i][sort_mac_idx[transcript][k + filled]] = r(k);
	  }
	  filled += m[transcript][j];
	}
	filled = 0;
	for (int j = mac_bins[transcript]; j < mac_bins[transcript] + maj_bins[transcript]; j++) {
	  arma::uvec r = unpack(tmp[j], m[transcript][j], j < mac_bins[transcript]); // Don't have to shuffle maj allele carriers
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
	build_minor_bins(odds, mac_indices, transcript, approximate, odds_sort, lower_bin_cutoff, upper_bin_cutoff);

	int msize = std::accumulate(m[transcript].begin(), m[transcript].end(), 0);
	assert(msize == mac_indices.n_elem);

	// Maj bin
	if(maj_indices.n_elem > 0) {
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
	  arma::uvec r = unpack(tmp[j], m[transcript][j], j < mac_bins[transcript]); // unpack and randomize cases into a vector
	  for (int k = 0; k < r.n_elem; k++) { //
		ret[transcript][i][sort_mac_idx[transcript][k + filled]] = r(k);
	  }
	  filled += m[transcript][j];
	}
	filled = 0;
	for (int j = mac_bins[transcript]; j < mac_bins[transcript] + maj_bins[transcript]; j++) {
	  arma::uvec r = unpack(tmp[j], m[transcript][j], j < mac_bins[transcript]); // Don't have to shuffle major allele carriers
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
  double nbins = approximate;
  double bin_width = ((arma::max(mac_odds_binnable) + 0.5) - arma::min(mac_odds_binnable)) / nbins;
  mac_bins[transcript] = nbins + lower.n_elem + upper.n_elem;
  double min_mac_odds = arma::min(mac_odds_binnable);
  // Fill in lower
  for(const auto &v : lower) {
	m[transcript].push_back(1);
	odds_[transcript].push_back(v);
  }
  for (arma::uword i = 0; i < nbins; i++) {
	std::vector<arma::uword> odds_in_range;
	for(arma::uword j = 0; j < mac_odds_binnable.n_elem; j++) { // sorted
	  if(mac_odds_binnable(j) >= min_mac_odds + i * bin_width && mac_odds_binnable(j) < min_mac_odds + (i + 1) * bin_width) {
		odds_in_range.push_back(j);
	  } else if(mac_odds_binnable(j) < min_mac_odds + i * bin_width) {
		continue;
	  } else {
		break;
	  }
	}
	arma::uvec odds_spanned = arma::conv_to<arma::uvec>::from(odds_in_range);
	if(odds_spanned.n_elem > 0) {
	  // Set number in group
	  m[transcript].push_back(odds_spanned.n_elem);
	  // Set odds for group
	  odds_[transcript].push_back(arma::mean(mac_odds_binnable(odds_spanned)));
	} else {
	  mac_bins[transcript]--;
	}
  }
  // Fill in upper
  for(const auto &v : upper) {
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
  double nbins = maj_nbins + lower.n_elem + upper.n_elem;
  maj_bins[transcript] = nbins;
  arma::uword stride = maj_odds_binnable.n_elem / nbins;
  arma::uvec odds_sort = arma::sort_index(maj_odds_binnable);
  sort_maj_idx[transcript] = maj_indices(arma::sort_index(maj_odds));
  // Fill in lower
  for(const auto &v : lower) {
	m[transcript].push_back(1);
	odds_[transcript].push_back(v);
  }
  for (arma::uword i = 0; i < maj_nbins; i++) {
	arma::span cur;
	if(i == maj_nbins - 1) { // Bin remaining, bin may be larger than other bins
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
  // Fill in upper
  for(const auto &v : upper) {
	m[transcript].push_back(1);
	odds_[transcript].push_back(v);
  }
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

	// Create bins
	sort_mac_idx[transcript] = arma::sort_index(odds); // Reuse sort_mac_idx for all bins
	arma::vec odds_sorted = arma::sort(odds);
	// subset bins
	double nbins = approximate;
	double bin_width = ((arma::max(odds) + 0.5) - arma::min(odds)) / nbins;
	mac_bins[transcript] = nbins;
	double min_odds = arma::min(odds);
	for (arma::uword i = 0; i < nbins; i++) {
	  std::vector<arma::uword> odds_in_range;
	  for(arma::uword j = 0; j < sort_mac_idx[transcript].n_elem; j++) { // sorted
		if(odds_sorted(j) >= min_odds + i * bin_width && odds_sorted(j) < min_odds + (i + 1) * bin_width) {
		  odds_in_range.push_back(j);
		} else if(odds_sorted(j) < min_odds + i * bin_width) {
		  continue;
		} else {
		  break;
		}
	  }
	  arma::uvec odds_spanned = arma::conv_to<arma::uvec>::from(odds_in_range);
	  if(odds_spanned.n_elem > 0) {
		// Set number in group
		m[transcript].push_back(odds_spanned.n_elem);
		// Set odds for group
		odds_[transcript].push_back(arma::mean(odds_sorted(odds_spanned)));
	  } else {
		mac_bins[transcript]--;
	  }
	}
	bins_built[transcript] = true;
  }

  for (int i = 0; i < nperm; i++) {
	std::vector<int32_t> tmp(odds_[transcript].size(), 0);
	sto.MultiFishersNCHyp(&tmp[0], &(m[transcript][0]), &(odds_[transcript][0]), ncases, odds_[transcript].size());

	// Unpack bins
	arma::uword filled = 0;
	for (int j = 0; j < mac_bins[transcript]; j++) {
	  arma::uvec r = unpack(tmp[j], m[transcript][j], true); // Always have to shuffle here
	  for (int k = 0; k < r.n_elem; k++) {
		  ret[transcript][i][sort_mac_idx[transcript](k + filled)] = r(k);
	  }
	  filled += m[transcript][j];
	}
  }
  return ret[transcript];
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
  for(arma::sword i = r.n_elem - 1; i >= 1; --i) {
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

std::vector<std::vector<int32_t>> Permute::permutations_mac_bin_fix(int nperm,
																	arma::vec &odds,
																	arma::uword ncases,
																	arma::uvec &mac_indices,
																	arma::uvec &maj_indices,
																	arma::uword &approximate,
																	const std::string &transcript) {
	if(bins_built.find(transcript) == bins_built.end()) {
		bins_built[transcript] = false;
	}

	arma::vec mac_odds = odds(mac_indices);
	arma::uvec mac_sort_idx = arma::sort_index(mac_odds);
	arma::uvec mac_unsort_idx = arma::sort_index(mac_sort_idx);
	arma::vec mac_sorted = mac_odds(mac_sort_idx);

	arma::vec maj_odds = odds(maj_indices);
	arma::uvec maj_sort_idx = arma::sort_index(maj_odds);
	arma::uvec maj_unsort_idx = arma::sort_index(maj_sort_idx);
	arma::vec maj_sorted = maj_odds(maj_sort_idx);

	if(!bins_built[transcript]) {
		ret[transcript] = std::vector<std::vector<int32_t>>(nperm);
		for (int i = 0; i < nperm; i++) {
			ret[transcript][i] = std::vector<int32_t>(odds.n_rows, 0);
		}

		// Create minor bins
		sort_mac_idx[transcript] = mac_indices(mac_sort_idx); // Not sure I need this anymore

		double nbins = std::max(approximate, 1ull);
		mac_bins[transcript] = nbins;
		double max_mac_sorted = arma::max(mac_sorted);
		double min_mac_sorted = arma::min(mac_sorted);
		double bin_width = ((max_mac_sorted + 0.5) - min_mac_sorted) / nbins;
		for (arma::uword i = 0; i < nbins; i++) {
			std::vector<arma::uword> odds_in_range;
			for(arma::uword j = 0; j < mac_sorted.n_elem; j++) { // sorted
				if(mac_sorted(j) >= min_mac_sorted + i * bin_width && mac_sorted(j) < min_mac_sorted + (i + 1) * bin_width) {
					odds_in_range.push_back(j);
				} else if(mac_sorted(j) < min_mac_sorted + i * bin_width) {
					continue;
				} else {
					break;
				}
			}
			arma::uvec odds_spanned = arma::conv_to<arma::uvec>::from(odds_in_range);
			if(odds_spanned.n_elem > 0) {
				// Set number in group
				m[transcript].push_back(odds_spanned.n_elem);
				// Set odds for group
				odds_[transcript].push_back(arma::mean(mac_sorted(odds_spanned)));
				mac_spans[transcript].push_back(odds_spanned);
			} else {
				mac_bins[transcript]--;
			}
		}

		// Maj bin
		// if(maj_indices.n_elem > 0) {
		// 	maj_bins[transcript] = 1;
		// 	sort_maj_idx[transcript] = maj_indices(maj_sort_idx); // Not sure I need this
		// 	// Set number in group
		// 	m[transcript].push_back(maj_indices.n_elem);
		// 	maj_spans[transcript] = arma::span::all;
		// 	// Set odds for group
		// 	odds_[transcript].push_back(arma::mean(odds(maj_indices)));
		// }
		if(maj_indices.n_elem > 0) {
			// nbins = std::min(approximate, maj_indices.n_elem);
            nbins = 10;
			maj_bins[transcript] = nbins;
			arma::uword stride = maj_indices.n_elem / nbins;
			arma::vec maj_odds = odds(maj_indices);
			arma::uvec odds_sort = arma::sort_index(maj_odds);
			sort_maj_idx[transcript] = maj_indices(odds_sort);
			for (arma::uword i = 0; i < nbins; i++) {
				arma::span cur(i * stride, std::min(i * stride + stride - 1, maj_indices.n_elem));

				maj_spans[transcript].push_back(cur);

				arma::vec odds_spanned = maj_odds(odds_sort(cur));
				// Set number in group
				m[transcript].push_back(odds_spanned.n_elem);
				// Set odds for group
				odds_[transcript].push_back(arma::mean(odds_spanned));
			}
		}
		bins_built[transcript] = true;
	}

	// Generate permutations
	for (int i = 0; i < nperm; i++) {
		arma::vec mac_fill(mac_indices.n_elem, arma::fill::zeros);
		arma::vec maj_fill(maj_indices.n_elem, arma::fill::zeros);
		std::vector<int32_t> tmp(odds_[transcript].size(), 0);
		sto.MultiFishersNCHyp(&tmp[0], &(m[transcript][0]), &(odds_[transcript][0]), ncases, odds_[transcript].size());

		// Unpack bins
		arma::uword filled = 0;
		for (int j = 0; j < mac_bins[transcript]; j++) { // for each bin
			arma::uvec r = unpack(tmp[j], m[transcript][j], j < mac_bins[transcript]); // unpack and randomize cases into a vector
			mac_fill(mac_spans[transcript][j]) = arma::conv_to<arma::vec>::from(r);
			// for (int k = 0; k < r.n_elem; k++) { //
			// 	ret[transcript][i][sort_mac_idx[transcript][k + filled]] = r(k);
			// }

			filled += m[transcript][j];
		}
		filled = 0;
		for (int j = mac_bins[transcript]; j < mac_bins[transcript] + maj_bins[transcript]; j++) {
			arma::uvec r = unpack(tmp[j], m[transcript][j], j < mac_bins[transcript]); // Don't have to shuffle major allele carriers
			maj_fill(maj_spans[transcript][j - mac_bins[transcript]]) = arma::conv_to<arma::vec>::from(r);
			// for (int k = 0; k < r.n_elem; k++) {
			// 	ret[transcript][i][sort_maj_idx[transcript][k + filled]] = r(k);
			// }
			filled += m[transcript][j];
		}
		arma::vec vtmp(odds.n_elem);
		vtmp(mac_indices) = mac_fill(mac_unsort_idx);
		vtmp(maj_indices) = maj_fill(maj_unsort_idx);

		for(int j = 0; j < ret[transcript][i].size(); j++) {
			ret[transcript][i][j] = vtmp(j);
		}
	}

	return ret[transcript];
};

