//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#include "carvaop.hpp"
#include "../utility/math.hpp"

CARVAOp::CARVAOp(CARVATask &ct, std::shared_ptr<Reporter> reporter, double seed, bool verbose)
	: done_(false),
	  verbose_(verbose),
	  carvaTask(ct),
	  reporter_(std::move(reporter)) {

}

CARVAOp::CARVAOp(CARVATask &&ct, std::shared_ptr<Reporter> reporter, double seed, bool verbose)
	: done_(false),
	  verbose_(verbose),
	  carvaTask(std::move(ct)),
	  reporter_(std::move(reporter)) {

}

auto CARVAOp::run() -> void {
  Stage stage = carvaTask.stage;

  if (stage == Stage::Stage1) {
	stage1();
  } else if (stage == Stage::Done) {
	carvaTask.methods.clear(carvaTask.gene.get_transcripts());
	this->done_ = true;
	return;
  } else {
	throw (std::runtime_error("Incorrect stage in CARVAOp::run()"));
  }
}

auto CARVAOp::finish() -> void {
  carvaTask.calc_multitranscript_pvalues();
  // Finalize result calculations and free memory
  carvaTask.cleanup();

  // Report results for non-gene_list run
  if (!carvaTask.tp.gene_list) {
	reporter_->sync_write_simple(carvaTask.results, carvaTask.tp);
	reporter_->sync_write_detail(carvaTask.gene.get_detail(), carvaTask.gene.is_testable());
	reporter_->sync_write_vaast(carvaTask, carvaTask.tp);
  }
}

auto CARVAOp::stage1() -> void {
  // Set original value
  for (auto &v : carvaTask.results) {
	std::string gene_name = v.first;
	std::string transcript = v.second.transcript;
	if (!carvaTask.gene.is_polymorphic(gene_name)) {
	  continue;
	}
	v.second.case_alt = arma::accu(carvaTask.gene.get_matrix(transcript).t() * carvaTask.get_cov().get_original_phenotypes());
	v.second.case_ref = 2 * arma::accu(carvaTask.get_cov().get_original_phenotypes()) - v.second.case_alt;
	v.second.cont_alt =
		arma::accu(carvaTask.gene.get_matrix(transcript).t() * (1. - carvaTask.get_cov().get_original_phenotypes()));
	v.second.cont_ref = 2 * arma::accu(1. - carvaTask.get_cov().get_original_phenotypes()) - v.second.case_alt;
	v.second.original =
		carvaTask.methods.call(carvaTask.gene, carvaTask.get_cov(), carvaTask.get_cov().get_original_phenotypes(), transcript, true);
  }

  int iter = 0;
  while (iter < carvaTask.npermutations) {
	int transcript_no = -1;

	for (auto &v : carvaTask.results) {
	  transcript_no++;
	  const std::string &k = v.second.transcript;
	  if (!carvaTask.gene.is_polymorphic(k)) {
		continue;
	  }

	  // Skip transcripts with no variants
	  if (!carvaTask.gene.is_polymorphic(k))
		continue;

	  arma::vec phenotypes = arma::conv_to<arma::vec>::from(carvaTask.get_permutations()[carvaTask.offset + iter]);
	  if (arma::accu(phenotypes) != arma::accu(carvaTask.cov->get_original_phenotypes())) {
		std::cerr << "vec count: " << std::accumulate(carvaTask.get_permutations()[carvaTask.offset + iter].begin(),
													  carvaTask.get_permutations()[carvaTask.offset + iter].end(),
													  0) << std::endl;
		std::cerr << "Permuted count: " << arma::accu(phenotypes) << std::endl;
		std::cerr << "Original count: " << arma::accu(carvaTask.cov->get_original_phenotypes()) << std::endl;
		throw (std::runtime_error("Failed to properly generate permutation."));
	  }

	  double perm_val = carvaTask.methods.call(carvaTask.gene, carvaTask.get_cov(), phenotypes, v.second.transcript, false);

	  v.second.permuted.push_back(perm_val);

	  // Update total number of permutations
	  v.second.permutations++;

	  check_perm(carvaTask.tp, perm_val, carvaTask.success_threshold, v);
	}
	// Stop iterating if everything is done
	if (std::all_of(carvaTask.results.cbegin(), carvaTask.results.cend(), [&](const auto &v) { return v.second.done; })) {
	  break;
	}
	iter++;
  }

  int ts_no = 0;
  for (auto &v : carvaTask.results) {
	double empirical;
	double midp;
	// If we stopped early, use the geometric correction, otherwise calculate for binomial phat
	if (carvaTask.tp.max_perms) {
	  if (v.second.permutations < *carvaTask.tp.max_perms) {
		empirical = geometric_p(v.second.successes, v.second.permutations);
		midp = geometric_p(v.second.mid_successes, static_cast<double>(v.second.permutations));
	  } else {
		empirical = (1. + v.second.successes) / (1. + v.second.permutations);
		midp = (1. + v.second.mid_successes) / (1. + v.second.permutations);
	  }
	} else {
	  if (v.second.permutations < carvaTask.tp.nperm) {
		empirical = geometric_p(v.second.successes, v.second.permutations);
		midp = geometric_p(v.second.mid_successes, static_cast<double>(v.second.permutations));
	  } else {
		empirical = (1. + v.second.successes) / (1. + v.second.permutations);
		midp = (1. + v.second.mid_successes) / (1. + v.second.permutations);
	  }
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
	if (carvaTask.tp.max_perms) {
	  if (carvaTask.tp.gene_list) {
		if (v.second.permutations >= carvaTask.termination) {
		  done_ = true;
		}
	  } else {
		if (v.second.permutations >= *carvaTask.tp.max_perms) {
		  done_ = true;
		}
	  }
	} else {
	  if (carvaTask.tp.gene_list) {
		if (v.second.permutations >= carvaTask.termination) {
		  done_ = true;
		}
	  } else {
		if (v.second.permutations >= carvaTask.tp.nperm) {
		  done_ = true;
		}
	  }
	}
	v.second.empirical_p = empirical;
	v.second.empirical_midp = midp;
	v.second.update_ci();

	double n = 2 * carvaTask.gene.get_samples().size();
	double nmac = arma::accu(arma::sum(arma::mat(carvaTask.gene.get_matrix(v.second.transcript)), 1) > 0);
	double nmaj = n - nmac;

	v.second.calc_exact_p(nmac, nmaj);
	ts_no++;
  }
  if (done_) { // True only when all are finished.
	if (verbose_) {
	  for (const auto &v : carvaTask.results) {
		std::cerr << carvaTask.gene.gene_name << "\t" << v.second.transcript << "\t";
		std::cerr << std::defaultfloat << std::setprecision(6) << v.second.original << std::endl;
	  }
	}
	carvaTask.stage = Stage::Done;
  } else {
	carvaTask.stage = Stage::Stage1;
  }
}

auto CARVAOp::check_perm(const TaskParams &tp,
						 double perm_val,
						 long success_threshold,
						 std::pair<const std::string, Result> &v) -> void {
  // Some methods return a pvalue so we need to reverse the success inequality
  if (tp.analytic) {
	if (perm_val <= v.second.original) {
	  if (tp.pthresh) {
		v.second.successes++;
		if (perm_val == v.second.original) {
		  v.second.mid_successes += 0.5;
		} else {
		  v.second.mid_successes++;
		}

		double lower, upper;
		std::tie(lower, upper) = poisson_ci(v.second.successes, v.second.permutations);

		// Ensure a minimum number of permutations
		if (v.second.permutations > 10) {
		  // Stop when the lower bound on the 95% confidence interval is greater than the threshold given
		  v.second.done = lower > *tp.pthresh;
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

		double lower, upper;
		std::tie(lower, upper) = poisson_ci(v.second.successes, v.second.permutations);

		// Ensure a minimum number of permutations
		if (v.second.permutations > 10) {
		  // Stop when the lower bound on the 95% confidence interval is greater than the threshold given
		  v.second.done = lower > *tp.pthresh;
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

