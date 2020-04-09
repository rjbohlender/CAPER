//
// Created by Bohlender,Ryan James on 2019-06-06.
//

#include "caeseop.hpp"
#include "../statistics/fishertest.hpp"

CAESEOp::CAESEOp(CAESETask &ct, std::shared_ptr<CAESEReporter> reporter, double seed, bool verbose)
	: gen_(seed), ta_(ct), done_(false), verbose_(verbose), reporter_(reporter) {

}

auto CAESEOp::run() -> void {
  effectsize();
}

auto CAESEOp::finish() -> void {
  ta_.cleanup();

  // Report results for non-gene_list run
  if (!ta_.get_tp().gene_list) {
	reporter_->sync_write_caese(ta_.results);
  }
}

auto CAESEOp::is_done() const -> bool {
  return done_;
}

auto CAESEOp::get_task() -> CAESETask {
  return ta_;
}

auto CAESEOp::effectsize() -> void {
  // Declare maps
  std::map<std::string, std::vector<int32_t>> mac_case_count;
  std::map<std::string, std::vector<int32_t>> maj_case_count;
  std::map<std::string, arma::colvec> mac_odds;
  std::map<std::string, arma::colvec> maj_odds;
  std::map<std::string, arma::uvec> mac_indices;
  std::map<std::string, arma::uvec> maj_indices;

  double perm_val;
  int transcript_no;
  // For permutation set output
  std::ofstream pset_ofs;


  // Setup
  for (auto &v : ta_.results) {
	const std::string &k = v.second.transcript;

	if (std::isnan(v.second.original)) {
	  FisherTest original(ta_.get_gene(), ta_.get_cov().get_original_phenotypes(), k);
	  v.second.original = original.get_or();
	  v.second.or_p = original.get_pval();
	  v.second.case_alt = original.case_alt;
	  v.second.case_ref = original.case_ref;
	  v.second.cont_alt = original.cont_alt;
	  v.second.cont_ref = original.cont_ref;
	}
	// Minor allele carrier indices
	mac_indices[k] = arma::find(arma::sum(arma::mat(ta_.get_gene().get_matrix(k)), 1) > 0);
	maj_indices[k] = arma::find(arma::sum(arma::mat(ta_.get_gene().get_matrix(k)), 1) == 0);
	assert(mac_indices[k].n_rows + maj_indices[k].n_rows == ta_.get_cov().get_nsamples());

	mac_odds[k] = ta_.get_cov().get_odds()(mac_indices[k]);
	maj_odds[k] = ta_.get_cov().get_odds()(maj_indices[k]);
  }

  int iter = 0;
  arma::vec phenotypes = ta_.get_cov().get_phenotype_vector();

  if (ta_.get_tp().permute_set) {
	pset_ofs.open(*ta_.get_tp().permute_set, std::ios_base::app);
  }

  while (iter < ta_.get_npermutations()) {
	// For each transcript in the gene
	transcript_no = -1;
	for (auto &v : ta_.results) {
	  std::vector<std::vector<int32_t>> permutations;
	  const std::string &k = v.second.transcript;

	  transcript_no++;

	  // Skip transcripts with no variants
	  if (!ta_.get_gene().is_polymorphic(k))
		continue;

	  // SKAT corrects for covariates so we don't use this permutation approach
	  if (ta_.get_tp().approximate) {
		permutations = ta_.get_permute(k).permutations_mac_bin(1,
															   ta_.get_cov().get_odds(),
															   ta_.get_cov().get_ncases(),
															   mac_indices[k],
															   maj_indices[k],
															   k,
															   *ta_.get_tp().approximate,
															   ta_.get_tp().maj_nbins,
															   ta_.get_tp().lower_bin_cutoff,
															   ta_.get_tp().upper_bin_cutoff);
	  } else {
		permutations = ta_.get_permute(k).permutations_maj_bin(1,
															   ta_.get_cov().get_odds(),
															   ta_.get_cov().get_ncases(),
															   mac_indices[k],
															   maj_indices[k],
															   k,
															   ta_.get_tp().maj_nbins,
															   ta_.get_tp().lower_bin_cutoff,
															   ta_.get_tp().upper_bin_cutoff);
	  }
	  phenotypes = arma::conv_to<arma::vec>::from(permutations[0]);

	  perm_val = call_method(ta_, phenotypes, ta_.get_tp(), k);

	  v.second.permuted.push_back(perm_val);

	  // Update total number of permutations
	  v.second.permutations++;
	  check_perm(ta_.get_tp(), perm_val, ta_.get_tp().success_threshold, v);

	  // Track when we reached threshold
	  if (v.second.successes == ta_.get_tp().success_threshold && v.second.min_success_at < 0) {
		v.second.min_success_at = v.second.permutations;
	  }
	}
	// Stop iterating if all transcripts are finished
	if (std::all_of(ta_.results.cbegin(), ta_.results.cend(), [&](const auto &v) { return v.second.done; }))
	  break;
	iter++;
  }
  if (verbose_) {
	if (!ta_.get_tp().alternate_permutation) {
	  for (const auto &v : ta_.results) {
	    double or_ = v.second.case_alt * v.second.cont_ref / (v.second.case_ref * v.second.cont_alt);
		std::cerr << "EffectSize: " << ta_.get_gene().get_gene() << "\t" << v.second.transcript << "\t";
		std::cerr << std::defaultfloat << std::setprecision(6) << v.second.original << "\t" << "vs. " << or_ << std::endl;
	  }
	} else {
	  for (const auto &v : ta_.results) {
		std::cerr << "Stage 2: " << ta_.get_gene().get_gene() << "\t" << v.second.transcript << "\t";
		std::cerr << std::defaultfloat << std::setprecision(6) << v.second.original << std::endl;
	  }
	}
  }
  if (ta_.get_tp().permute_set) {
	pset_ofs.close();
	std::exit(0);
  }

  for (auto &v : ta_.results) {
	double empirical;
	double midp;
	if (v.second.min_success_at > 0 && v.second.successes == ta_.get_tp().success_threshold + 1) {
	  std::uniform_int_distribution<> dis(v.second.min_success_at, v.second.permutations);

	  empirical = v.second.successes / (1. + dis(gen_));
	  midp = v.second.mid_successes / (1 + dis(gen_));
	} else {
	  empirical = (1. + v.second.successes) / (1. + v.second.permutations);
	  midp = (1. + v.second.mid_successes) / (1. + v.second.permutations);
	}

	// Success on every iteration
	if (empirical > 1) {
	  empirical = 1;
	}

	if (midp > 1) {
	  midp = 1;
	}

	v.second.empirical_p = empirical;
	v.second.empirical_midp = midp;
  }
  done_ = true;
}

auto CAESEOp::check_perm(const TaskParams &tp,
						 double perm_val,
						 int success_threshold,
						 std::pair<const std::string, Result> &v) -> void {
  if (perm_val >= v.second.original) {
	if (tp.pthresh) {
	  v.second.successes++;
	  if (perm_val == v.second.original) {
		v.second.mid_successes += 0.5;
	  } else {
		v.second.mid_successes++;
	  }

	  double p = static_cast<double>(v.second.successes) / static_cast<double>(v.second.permutations);

	  // 95% confidence interval
	  double ci = 1.96 * std::sqrt(p * (1. - p) / v.second.permutations);
	  // Ensure a minimum number of permutations
	  if (v.second.permutations > 10) {
		// Stop when the lower bound on the 95% confidence interval is greater than the threshold given
		v.second.done = p - ci > *tp.pthresh;
	  }
	} else if (v.second.successes < success_threshold) {
	  v.second.successes++;
	  if (perm_val == v.second.original) {
		v.second.mid_successes += 0.5;
	  } else {
		v.second.mid_successes++;
	  }
	} else {
	  v.second.successes++;
	  if (perm_val == v.second.original) {
		v.second.mid_successes += 0.5;
	  } else {
		v.second.mid_successes++;
	  }
	  v.second.done = true;
	}
  }
}

auto CAESEOp::call_method(CAESETask &ct, arma::vec &phenotypes, TaskParams &tp, const std::string &k) -> double {
  return FisherTest(ct.get_gene(), phenotypes, k).get_or();
}
