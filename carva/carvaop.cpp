//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#include "carvaop.hpp"

CARVAOp::CARVAOp(CARVATask &ta, std::shared_ptr<Reporter> reporter, double seed, bool verbose)
: ta_(ta),
  gen_(seed),
  reporter_(reporter),
  done_(false),
  verbose_(verbose) {

}

CARVAOp::CARVAOp(const CARVAOp &op)
: ta_(op.ta_),
  done_(op.done_),
  verbose_(op.verbose_),
  reporter_(op.reporter_) {
}

CARVAOp::CARVAOp(CARVAOp &&op) noexcept
: ta_(op.ta_),
  done_(op.done_),
  verbose_(op.verbose_),
  reporter_(op.reporter_) {

}

CARVAOp &CARVAOp::operator=(const CARVAOp &rhs) {
  ta_ = rhs.ta_;
  done_ = rhs.done_;
  verbose_ = rhs.verbose_;
  reporter_ = rhs.reporter_;

  return *this;
}

auto CARVAOp::run() -> void {
  Stage stage = ta_.get_stage();

  if (stage == Stage::Stage1) {
	stage1();
  } else if (stage == Stage::Stage2) {
	stage2();
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
  }
}

auto CARVAOp::stage1() -> void {
  // Set original value
  for (auto &v : ta_.results) {
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
						false,
						false);
	  } else {
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

	  // Update tota_l number of permuta_tions
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
	}
  }
}

auto CARVAOp::stage2() -> void {
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
  if (!ta_.get_tp().alternate_permutation) {
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
	  // Minor allele carrier indices
	  mac_indices[k] = arma::find(arma::sum(arma::mat(ta_.get_gene().get_matrix(k)), 1) > 0);
	  maj_indices[k] = arma::find(arma::sum(arma::mat(ta_.get_gene().get_matrix(k)), 1) == 0);
	  assert(mac_indices[k].n_rows + maj_indices[k].n_rows == ta_.get_cov().get_nsamples());

	  mac_odds[k] = ta_.get_cov().get_odds()(mac_indices[k]);
	  maj_odds[k] = ta_.get_cov().get_odds()(maj_indices[k]);
	}
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

	  transcript_no++;

	  // Skip transcripts with no variants
	  if (!ta_.get_gene().is_polymorphic(k))
		continue;

	  // SKAT corrects for covariates so we don't use this permutation approach
	  if (!ta_.get_tp().alternate_permutation) {
		if(ta_.get_tp().linear) {
		  // Fisher-Yates Shuffle
		  for(arma::sword i = phenotypes.n_elem - 1; i > 0; --i) {
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
																  mac_indices[k],
																  maj_indices[k],
																  *ta_.get_tp().approximate,
																  k);
		  } else {
			permutations = ta_.get_permute(k).permutations_maj_bin(1,
																  ta_.get_cov().get_odds(),
																  ta_.get_cov().get_ncases(),
																  mac_indices[k],
																  maj_indices[k],
																  k);
		  }
		  phenotypes = arma::conv_to<arma::vec>::from(permutations[0]);
		}

		perm_val = call_method(ta_.get_methods(),
							   ta_.get_gene(),
							   ta_.get_cov(),
							   phenotypes,
							   ta_.get_tp(),
							   k,
							   false,
							   false);
	  } else {
		perm_val = call_method(ta_.get_methods(),
							   ta_.get_gene(),
							   ta_.get_cov(),
							   phenotypes,
							   ta_.get_tp(),
							   k,
							   transcript_no == 0,
							   false);
	  }

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
	if (!ta_.get_tp().alternate_permutation) {
	  for (const auto &v : ta_.results) {
		std::cerr << "Stage 2: " << ta_.get_gene().get_gene() << "\t" << v.second.transcript << "\t";
		std::cerr << std::defaultfloat << std::setprecision(6) << v.second.original << std::endl;
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
	if (v.second.min_success_at > 0 && v.second.successes == ta_.success_threshold + 1) {
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
  ta_.set_stage(Stage::Done);
  done_ = true;
}

auto CARVAOp::is_done() const -> bool {
  return done_;
}

auto CARVAOp::check_perm(const TaskParams &tp,
						  double perm_val,
						  int success_threshold,
						  std::pair<const std::string, Result> &v) -> void {
  // SKATO returns a pvalue so we need to reverse the successes
  if (tp.method == "SKATO" || (tp.method == "SKAT" && tp.total_permutations == 0) || tp.method == "CMC" || tp.method == "RVT1" || tp.method == "RVT2") {
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
	return method.BURDEN(gene, k, shuffle, tp.a, tp.b);
  } else if (tp.method == "CALPHA") {
	return method.CALPHA(gene, phenotypes, k);
  } else if (tp.method == "CMC") {
	return method.CMC(gene, phenotypes, k, tp.cmcmaf);
  } else if (tp.method == "RVT1") {
	return method.RVT1(gene, phenotypes, cov.get_covariate_matrix(), k, tp.linear);
  } else if (tp.method == "RVT2") {
	return method.RVT2(gene, phenotypes, cov.get_covariate_matrix(), k, tp.linear);
  } else if (tp.method == "SKAT") {
	return method.SKATR(gene, k, shuffle, tp.a, tp.b, detail, tp.linear, tp.total_permutations > 0);
  } else if (tp.method == "SKATO") {
	return method.SKATRO(gene, k, shuffle, tp.a, tp.b, detail, tp.linear);
  } else if (tp.method == "VAAST") {
	return method.Vaast(gene,
						phenotypes,
						k,
						tp.score_only_minor,
						tp.score_only_alternative,
						2.0,
						tp.group_size,
						detail,
						tp.biallelic);
  } else if (tp.method == "VT") {
	return method.VT(gene, k, shuffle);
  } else if (tp.method == "WSS") {
	return method.WSS(gene, phenotypes, k);
  } else {
	throw(std::logic_error("Failed to find method."));
  }
}

auto CARVAOp::get_args() -> CARVATask {
  return ta_;
}

