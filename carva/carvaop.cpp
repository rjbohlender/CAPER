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
  Stage stage = ta_.stage;

  if (stage == Stage::Stage1) {
	stage1();
  } else if (stage == Stage::Done) {
	ta_.methods.clear(ta_.gene.get_transcripts());
	this->done_ = true;
	return;
  } else {
	throw (std::runtime_error("Incorrect stage in CARVAOp::run()"));
  }
}

auto CARVAOp::finish() -> void {
  ta_.calc_multitranscript_pvalues();
  // Finalize result calculations and free memory
  ta_.cleanup();

  // Report results for non-gene_list run
  if (!ta_.tp.gene_list) {
	reporter_->sync_write_simple(ta_.results, ta_.tp);
	reporter_->sync_write_detail(ta_.gene.get_detail(), ta_.gene.is_testable());
	reporter_->sync_write_vaast(ta_, ta_.tp);
  }
}

auto CARVAOp::stage1() -> void {
  // Set original value
  for (auto &v : ta_.results) {
    std::string gene_name = v.first;
    std::string transcript = v.second.transcript;
	if (!ta_.gene.is_polymorphic(gene_name)) {
	  continue;
	}
	v.second.case_alt = arma::accu(ta_.gene.get_matrix(transcript).t() * ta_.get_cov().get_original_phenotypes());
	v.second.case_ref = 2 * arma::accu(ta_.get_cov().get_original_phenotypes()) - v.second.case_alt;
	v.second.cont_alt = arma::accu(ta_.gene.get_matrix(transcript).t() * (1. - ta_.get_cov().get_original_phenotypes()));
	v.second.cont_ref = 2 * arma::accu(1. - ta_.get_cov().get_original_phenotypes()) - v.second.case_alt;
	v.second.original =
		call_method(ta_.methods,
					ta_.gene,
					ta_.get_cov(),
					ta_.get_cov().get_original_phenotypes(),
					ta_.tp,
					transcript,
					false,
					true);
  }

  int iter = 0;
  while (iter < ta_.npermutations) {
	int transcript_no = -1;

	for (auto &v : ta_.results) {
	  transcript_no++;
	  const std::string &k = v.second.transcript;
	  if (!ta_.gene.is_polymorphic(k)) {
		continue;
	  }

	  double perm_val;

	  // Skip transcripts with no variants
	  if (!ta_.gene.is_polymorphic(k))
		continue;

	  arma::vec phenotypes;
	  phenotypes = arma::conv_to<arma::vec>::from(ta_.get_permutations()[iter]);
	  perm_val = call_method(ta_.methods,
							 ta_.gene,
							 ta_.get_cov(),
							 phenotypes,
							 ta_.tp,
							 v.second.transcript,
							 transcript_no == 0,
							 false);

	  // ta_.increment_permuted(v.second.transcript, perm_val);
	  v.second.permuted.push_back(perm_val);

	  // Update total number of permutations
	  v.second.permutations++;

	  check_perm(ta_.tp, perm_val, ta_.success_threshold, v);

	  // Track when we reached success threshold
	  if (v.second.successes == ta_.success_threshold && v.second.min_success_at < 0) {
		v.second.min_success_at = v.second.permutations;
	  }
	}
	// Stop iterating if everything is done
	if (std::all_of(ta_.results.cbegin(), ta_.results.cend(), [&](const auto &v) { return v.second.done; })) {
	  break;
	}
	iter++;
  }

  int ts_no = 0;
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

	if (ts_no == 0) {
	  done_ = v.second.done;
	} else {
	  done_ &= v.second.done;
	}
	if(ta_.tp.max_perms) {
	  if (ta_.tp.gene_list) {
		if (v.second.permutations >= *ta_.tp.max_perms / (ta_.tp.nthreads - 1)) {
		  done_ = true;
		}
	  } else {
		if (v.second.permutations >= *ta_.tp.max_perms) {
		  done_ = true;
		}
	  }
	} else {
	  if (ta_.tp.gene_list) {
		if (v.second.permutations >= ta_.tp.nperm / (ta_.tp.nthreads - 1)) {
		  done_ = true;
		}
	  } else {
		if (v.second.permutations >= ta_.tp.nperm) {
		  done_ = true;
		}
	  }
	}
	v.second.empirical_p = empirical;
	v.second.empirical_midp = midp;
	v.second.update_ci();

	double n = 2 * ta_.gene.get_samples().size();
	double nmac = arma::accu(arma::sum(arma::mat(ta_.gene.get_matrix(v.second.transcript)), 1) > 0);
	double nmaj = n - nmac;

	v.second.calc_exact_p(nmac, nmaj);
	ts_no++;
  }
  if (done_) { // True only when all are finished.
	if (verbose_) {
	  for (const auto &v : ta_.results) {
		std::cerr << ta_.gene.gene_name << "\t" << v.second.transcript << "\t";
		std::cerr << std::defaultfloat << std::setprecision(6) << v.second.original << std::endl;
	  }
	}
    ta_.stage = Stage::Done;
  } else {
    ta_.stage = Stage::Stage1;
  }
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
						  const TaskParams &tp,
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
	return method.SKAT(gene, k, phenotypes, tp.a, tp.b, detail, tp.linear, tp.nperm > 0, true);
  } else if (tp.method == "SKATO") {
	return method.SKATO(gene, k, phenotypes, tp.a, tp.b, detail, tp.linear);
  } else if (tp.method == "VAAST") {
	return method.VAAST(gene,
						phenotypes,
						k,
						tp.vaast_site_penalty,
						tp.group_size,
						detail,
						tp.biallelic,
						tp.soft_maf_filter,
						tp.legacy_grouping);
  } else if (tp.method == "VT") {
	return method.VT(gene, k, phenotypes);
  } else if (tp.method == "WSS") {
	return method.WSS(gene, phenotypes, k);
  } else {
	throw (std::logic_error("Failed to find method."));
  }
}
