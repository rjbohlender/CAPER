//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#include "carvaop.hpp"

CARVAOp::CARVAOp(CARVATask &ct, std::shared_ptr<Reporter> reporter, double seed, bool verbose)
	: gen_(seed),
	  ta_(ct),
	  done_(false),
	  verbose_(verbose),
	  reporter_(std::move(reporter)) {

}

CARVAOp::CARVAOp(CARVATask &&ct, std::shared_ptr<Reporter> reporter, double seed, bool verbose)
	: gen_(seed),
	  ta_(std::move(ct)),
	  done_(false),
	  verbose_(verbose),
	  reporter_(std::move(reporter)) {

}

auto CARVAOp::run() -> void {
  Stage stage = ta_.get_stage();

  if (stage == Stage::Stage1) {
	stage1();
  } else if (stage == Stage::Stage2) {
	stage2();
  } else if (stage == Stage::Done) {
    ta_.get_methods().clear(ta_.get_gene().get_transcripts());
    ta_.get_cov().clear();
    return;
  } else {
    throw(std::runtime_error("Incorrect stage in CARVAOp::run()"));
  }
}

auto CARVAOp::finish() -> void {
  ta_.calc_multitranscript_pvalues();
  // Finalize result calculations and free memory
  ta_.cleanup();

  // Report results for non-gene_list run
  if (!ta_.get_tp().gene_list) {
    reporter_->sync_write_simple(ta_.results, ta_.get_tp(), ta_.get_tp().top_only);
    reporter_->sync_write_detail(ta_.get_gene().get_detail(), ta_.get_gene().is_testable());
    reporter_->sync_write_vaast(ta_, ta_.get_tp());
  }
}

auto CARVAOp::stage1() -> void {
  // Set original value
  for (auto &v : ta_.results) {
    if(!ta_.get_gene().is_polymorphic(v.first)) {
      continue;
    }
	v.second.original =
		call_method(ta_.get_methods(),
					ta_.get_gene(),
					ta_.get_cov(),
					ta_.get_cov().get_original_phenotypes(),
					ta_.get_tp(),
					v.second.transcript,
					false,
					true);
  }


  if (verbose_) {
	if (!ta_.get_tp().alternate_permutation) {
	  for (const auto &v : ta_.results) {
		std::cerr << "Stage 1: " << ta_.get_gene().get_gene() << "\t" << v.second.transcript << "\t";
		std::cerr << std::defaultfloat << std::setprecision(6) << v.second.original << std::endl;
	  }
	} else {
	  for (const auto &v : ta_.results) {
		std::cerr << "Stage 1: " << ta_.get_gene().get_gene() << "\t" << v.second.transcript << "\t";
		std::cerr << std::defaultfloat << std::setprecision(6) << v.second.original << std::endl;
	  }
	}
  }

  int iter = 0;
  while (iter < ta_.get_npermutations()) {
	int transcript_no = -1;

	for (auto &v : ta_.results) {
	  transcript_no++;
	  const std::string &k = v.second.transcript;
      if(!ta_.get_gene().is_polymorphic(k)) {
        continue;
      }

	  double perm_val;

	  // Skip transcripts with no variants
	  if (!ta_.get_gene().is_polymorphic(k))
		continue;

      arma::vec phenotypes;
	  if (!ta_.get_tp().alternate_permutation) {
		phenotypes = arma::conv_to<arma::vec>::from(ta_.get_permutations()[iter]);
		perm_val =
			call_method(ta_.get_methods(),
						ta_.get_gene(),
						ta_.get_cov(),
						phenotypes,
						ta_.get_tp(),
						v.second.transcript,
						transcript_no == 0,
						false);
	  } else {
		phenotypes = arma::conv_to<arma::vec>::from(ta_.get_permutations()[iter]);
		perm_val = call_method(ta_.get_methods(),
							   ta_.get_gene(),
							   ta_.get_cov(),
							   phenotypes,
							   ta_.get_tp(),
							   v.second.transcript,
							   transcript_no == 0,
							   false);
	  }

	  // ta_.increment_permuted(v.second.transcript, perm_val);
	  v.second.permuted.push_back(perm_val);

	  // Update total number of permutations
	  v.second.permutations++;

	  check_perm(ta_.get_tp(), perm_val, ta_.success_threshold, v);

	  // Track when we reached success threshold
	  if (v.second.successes == ta_.success_threshold && v.second.min_success_at < 0) {
		v.second.min_success_at = v.second.permutations;
	  }
	}
	// Stop iterating if everything is done
	if (std::all_of(ta_.results.cbegin(), ta_.results.cend(), [&](const auto &v) { return v.second.done; }))
	  break;
	iter++;
  }

  if (std::any_of(ta_.results.cbegin(), ta_.results.cend(), [&](const auto &v) { return !v.second.done; })
	  && ta_.get_remaining() > 0) {
	ta_.set_stage(Stage::Stage2);
  } else {
	for (auto &v : ta_.results) {
	  double empirical;
	  double midp;
	  if (v.second.min_success_at > 0 && v.second.successes == ta_.success_threshold + 1) {
		std::uniform_int_distribution<> dis(v.second.min_success_at, v.second.permutations);

		empirical = v.second.successes / (1. + dis(gen_));
		midp = v.second.mid_successes / (1. + dis(gen_));
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

	  ta_.set_stage(Stage::Done);
	  done_ = true;
	  v.second.empirical_p = empirical;
	  v.second.empirical_midp = midp;
	  v.second.update_ci();

	  double n = 2 * ta_.get_gene().get_samples().size();
	  double nmac = arma::accu(arma::sum(arma::mat(ta_.get_gene().get_matrix(v.second.transcript)), 1) > 0);
	  double nmaj = n - nmac;

	  v.second.calc_exact_p(nmac, nmaj);
	}
  }
}

auto CARVAOp::stage2() -> void {
  arma::vec mac_odds;
  arma::vec maj_odds;
  arma::uvec mac_indices;
  arma::uvec maj_indices;

  double perm_val = 0;
  int transcript_no = -1;
  // For permutation set output
  std::ofstream pset_ofs;

  // Setup
  if (!ta_.get_tp().alternate_permutation) {
    arma::vec mac_carriers(ta_.get_cov().get_nsamples(), arma::fill::zeros);
	for (auto &v : ta_.results) {
	  const std::string &k = v.second.transcript;
	  if(!ta_.get_gene().is_polymorphic(k)) {
	    continue;
	  }

	  if (std::isnan(v.second.original)) {
		v.second.original = call_method(ta_.get_methods(),
										ta_.get_gene(),
										ta_.get_cov(),
										ta_.get_cov().get_original_phenotypes(),
										ta_.get_tp(),
										k,
										false,
										true);
	  }
	  mac_carriers(arma::find(arma::sum(arma::mat(ta_.get_gene().get_matrix(k) + ta_.get_gene().get_missing(k)), 1) > 0)).ones();
	}
	// Minor and major allele carrier indices
	mac_indices = arma::find(mac_carriers > 0);
	// Missing are permuted with mac
	maj_indices = arma::find(mac_carriers == 0);
	// Missing are excluded from maj
	assert(mac_indices.n_rows + maj_indices.n_rows == ta_.get_cov().get_nsamples());

	mac_odds = ta_.get_cov().get_odds()(mac_indices);
	maj_odds = ta_.get_cov().get_odds()(maj_indices);
  } else {
	for (auto &v : ta_.results) {
	  const std::string &k = v.second.transcript;

	  if (std::isnan(v.second.original)) {
		v.second.original = call_method(ta_.get_methods(),
										ta_.get_gene(),
										ta_.get_cov(),
										ta_.get_cov().get_original_phenotypes(),
										ta_.get_tp(),
										 k,
										false,
										true);
	  }
	}
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

      if(!ta_.get_gene().is_polymorphic(k)) {
        continue;
      }

	  transcript_no++;

	  // Skip transcripts with no variants
	  if (!ta_.get_gene().is_polymorphic(k))
		continue;

	  if (transcript_no == 0) { // Run permutation for only the first transcript and reuse to maintain MGIT functionality
		if (ta_.get_tp().linear || !ta_.get_tp().covadj) {
		  // Fisher-Yates Shuffle
		  for (arma::sword i = phenotypes.n_elem - 1; i > 0; --i) {
			std::uniform_int_distribution<> dis(0, i);
			auto j = static_cast<arma::uword>(dis(gen_));
			double tmp = phenotypes(i);
			phenotypes(i) = phenotypes(j);
			phenotypes(j) = tmp;
		  }
		} else {
		  if (ta_.get_tp().approximate) {
			permutations = ta_.get_permute(k).permutations_mac_bin(1,
																   ta_.get_cov().get_odds(),
																   ta_.get_cov().get_ncases(),
																   mac_indices,
																   maj_indices,
																   k,
																   *ta_.get_tp().approximate,
																   ta_.get_tp().maj_nbins,
																   ta_.get_tp().lower_bin_cutoff,
																   ta_.get_tp().upper_bin_cutoff);
		  } else {
			permutations = ta_.get_permute(k).permutations_maj_bin(1,
																   ta_.get_cov().get_odds(),
																   ta_.get_cov().get_ncases(),
																   mac_indices,
																   maj_indices,
																   k,
																   ta_.get_tp().maj_nbins,
																   ta_.get_tp().lower_bin_cutoff,
																   ta_.get_tp().upper_bin_cutoff);
		  }
		  phenotypes = arma::conv_to<arma::vec>::from(permutations[0]);
		  if(ta_.get_tp().permute_set) {
			pset_ofs << phenotypes.t();
		  }
		}
	  }
	  perm_val = call_method(ta_.get_methods(),
							 ta_.get_gene(),
							 ta_.get_cov(),
							 phenotypes,
							 ta_.get_tp(),
							  k,
							 true,
							 false);

	  // ta.increment_permuted(v.second.transcript, perm_val);
	  v.second.permuted.push_back(perm_val);

	  // Update total number of permutations
	  v.second.permutations++;
	  check_perm(ta_.get_tp(), perm_val, ta_.success_threshold, v);

	  // Track when we reached threshold

	  if (v.second.successes == ta_.success_threshold && v.second.min_success_at < 0) {
		v.second.min_success_at = v.second.permutations;
	  }
	}
	// Stop iterating if all transcripts are finished
	if (std::all_of(ta_.results.cbegin(), ta_.results.cend(), [&](const auto &v) { return v.second.done; }))
	  break;
	iter++;
  }
  if (verbose_) {
	for (const auto &v : ta_.results) {
	  std::cerr << "Stage 2: " << ta_.get_gene().get_gene() << "\t" << v.second.transcript << "\t";
	  std::cerr << std::defaultfloat << std::setprecision(6) << v.second.original << std::endl;
	}
  }
  if (ta_.get_tp().permute_set) {
	pset_ofs.close();
	std::exit(0);
  }

  for (auto &v : ta_.results) {
	double empirical;
	double midp;
	if (v.second.min_success_at > 0 && v.second.successes == ta_.success_threshold + 1) {
	  std::uniform_int_distribution<> dis(v.second.min_success_at, v.second.permutations);

	  int rand_perms = dis(gen_);
	  v.second.rand_perms = rand_perms;

	  empirical = v.second.successes / (1. + rand_perms);
	  midp = v.second.mid_successes / (1 + rand_perms);
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
	v.second.update_ci();

	double n = 2 * ta_.get_gene().get_samples().size();
	double nmac = arma::accu(arma::sum(arma::mat(ta_.get_gene().get_matrix(v.second.transcript)), 1) > 0);
	double nmaj = n - nmac;

	v.second.calc_exact_p(nmac, nmaj);
  }
  ta_.set_stage(Stage::Done);
  done_ = true;
}

auto CARVAOp::check_perm(const TaskParams &tp,
						 double perm_val,
						 int success_threshold,
						 std::pair<const std::string, Result> &v) -> void {
  // SKATO returns a pvalue so we need to reverse the successes
  if (tp.analytic) {
	if (perm_val <= v.second.original) {
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
  } else {
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
}

auto CARVAOp::call_method(Methods &method,
						  Gene &gene,
						  Covariates &cov,
						  arma::vec &phenotypes,
						  TaskParams &tp,
						  const std::string &k,
						  bool shuffle,
						  bool detail) -> double {
  if (tp.method == "BURDEN") {
	return method.BURDEN(gene, k, phenotypes, tp.a, tp.b);
  } else if (tp.method == "CALPHA") {
	return method.CALPHA(gene, phenotypes, k);
  } else if (tp.method == "CMC") {
	return method.CMC(gene, phenotypes, k, tp.cmcmaf);
  } else if (tp.method == "CMC1df") {
    return method.CMC1df(gene, phenotypes, k);
  } else if (tp.method == "RVT1") {
	return method.RVT1(gene, phenotypes, cov.get_covariate_matrix(), cov.get_coef(), k, tp.linear);
  } else if (tp.method == "RVT2") {
	return method.RVT2(gene, phenotypes, cov.get_covariate_matrix(), cov.get_coef(), k, tp.linear);
  } else if (tp.method == "SKAT") {
	return method.SKAT(gene, k, phenotypes, tp.a, tp.b, detail, tp.linear, tp.total_permutations > 0, true);
  } else if (tp.method == "SKATO") {
	return method.SKATO(gene, k, phenotypes, tp.a, tp.b, detail, tp.linear);
  } else if (tp.method == "VAAST") {
	return method.VAAST(gene,
						phenotypes,
						k,
						tp.score_only_minor,
						tp.score_only_alternative,
						tp.vaast_site_penalty,
						tp.group_size,
						detail,
						tp.biallelic,
						tp.soft_maf_filter);
  } else if (tp.method == "VT") {
	return method.VT(gene, k, phenotypes);
  } else if (tp.method == "WSS") {
	return method.WSS(gene, phenotypes, k);
  } else {
	throw (std::logic_error("Failed to find method."));
  }
}
