//
// Created by Bohlender,Ryan James on 2019-03-19.
//

#include "powerop.hpp"
#include <iostream>
#include <iomanip>

PowerOp::PowerOp(PowerTask &pt, std::shared_ptr<PowerReporter> reporter, double seed, bool verbose)
	: pt_(pt), reporter_(reporter), gen_(seed), bootstrapped_(pt_.gene), done_(false) {
}

auto PowerOp::run() -> void {
  power();
}

auto PowerOp::finish() -> void {
  // cleanup
  if (!pt_.tp.gene_list) {
    reporter_->sync_write_power(pt_.power_res_);
  }
}

auto PowerOp::is_done() const -> bool {
  return done_;
}

auto PowerOp::get_task() -> PowerTask {
  return pt_;
}

auto PowerOp::power() -> void {
  if (cases_.empty()) {
	cases_ = arma::find(pt_.cov->get_original_phenotypes() > 0);
  }
  if (controls_.empty()) {
	controls_ = arma::find(pt_.cov->get_original_phenotypes() == 0);
  }

  for (const auto &ts : pt_.gene.get_transcripts()) {
	if (pt_.tp.alpha.size() > 1) {
	  for (arma::uword k = 0; k < pt_.tp.ncases.size(); k++) {
		// Alpha vector to track the successes
		success_map_[ts].push_back(arma::vec(pt_.tp.alpha.size(), arma::fill::zeros));
		arma::uword ncases = pt_.tp.ncases[k];
		arma::uword ncontrols = pt_.tp.ncontrols[k];

		phenotypes_.zeros(ncases + ncontrols);
		phenotypes_(arma::span(0, ncases - 1)).fill(1);

		arma::vec odds = pt_.cov->get_odds();
		for (arma::uword i = 0; i < pt_.nreps; i++) {
		  bootstrapped_.set_matrix(ts, sample(pt_.gene.get_matrix(ts), ncases, ncontrols));
		  // Original for this bootstrap replicate
		  double original = call_method(pt_.method, bootstrapped_, *pt_.cov, phenotypes_, pt_.tp, ts);
		  // Need to get a p-value via permutation for some methods
		  // If the method returns a p-value
		  if (pt_.tp.method == "SKATO" || pt_.tp.method == "SKAT" || pt_.tp.method == "CMC" || pt_.tp.method == "RVT1"
			  || pt_.tp.method == "RVT2") {
			double val = call_method(pt_.method, bootstrapped_, *pt_.cov, phenotypes_, pt_.tp, ts);
			if (val < original) {
			  success_map_[ts][k]++;
			}
		  } else {
			// We need to permute to get a p-value
			// Cease permutation if the p-value ci excludes alpha above or below
			// If below, call it a success
			int block = static_cast<int>(pt_.tp.nperm); // Permutation block size
			double successes = 0;
			double p;
			double val;

			arma::uvec idx = arma::join_vert(cases_(case_idx_), controls_(control_idx_));
			arma::vec bootstrapped_odds = odds(idx);

			// Shuffle phenotypes
			std::vector<std::vector<int8_t>> permutations = permute_.epsilon_permutation(
				block,
				bootstrapped_odds,
				ncases,
				ts,
				pt_.tp.bin_epsilon);

			for (arma::uword j = 0; j < block; j++) {
			  phenotypes_ = arma::conv_to<arma::vec>::from(permutations[j]);
			  pt_.cov->set_phenotype_vector(phenotypes_);
			  val = call_method(pt_.method, bootstrapped_, *pt_.cov, phenotypes_, pt_.tp, ts);
			  if (val >= original) {
				successes++;
			  }
			}

			p = (successes + 2.) / (block + 4.); // Wilson Score estimate
			success_map_[ts][k](arma::find(pt_.tp.alpha >= p)) += 1;
		  }
		}
	  }
	} else {
	  // Single alpha value,
	  for (const auto &a : pt_.tp.alpha) {
		for (arma::uword k = 0; k < pt_.tp.ncases.size(); k++) {
		  arma::uword ncases = pt_.tp.ncases[k];
		  arma::uword ncontrols = pt_.tp.ncontrols[k];

		  phenotypes_.zeros(ncases + ncontrols);
		  phenotypes_(arma::span(0, ncases - 1)).fill(1);

		  arma::vec odds = pt_.cov->get_odds();
		  for (arma::uword i = 0; i < pt_.nreps; i++) {
			bootstrapped_.set_matrix(ts, sample(pt_.gene.get_matrix(ts), ncases, ncontrols));
			// Original for this bootstrap replicate
			double original = call_method(pt_.method, bootstrapped_, *pt_.cov, phenotypes_, pt_.tp, ts);
			// Need to get a p-value via permutation for some methods
			// If the method returns a p-value
			if (pt_.tp.method == "SKATO" || pt_.tp.method == "SKAT" || pt_.tp.method == "CMC" || pt_.tp.method == "RVT1"
				|| pt_.tp.method == "RVT2") {
			  double val = call_method(pt_.method, bootstrapped_, *pt_.cov, phenotypes_, pt_.tp, ts);
			  if (val < original) {
				success_map_[ts][k]++;
			  }
			} else {
			  // We need to permute to get a p-value
			  // Cease permutation if the p-value ci excludes alpha above or below
			  // If below, call it a success
			  arma::uword n = 0;
			  int block = 5; // Permutation block size
			  double successes = 0;
			  double p, sp;
			  double z, z2;
			  double val;
			  double ci;

			  arma::uvec idx = arma::join_vert(cases_(case_idx_), controls_(control_idx_));
			  arma::vec bootstrapped_odds = odds(idx);
			  do {
				// Shuffle phenotypes
				std::vector<std::vector<int8_t>> permutations = permute_.epsilon_permutation(
					block,
					bootstrapped_odds,
					ncases,
					ts,
					pt_.tp.bin_epsilon);

				for (arma::uword j = 0; j < block; j++) {
				  n++;
				  phenotypes_ = arma::conv_to<arma::vec>::from(permutations[j]);
				  pt_.cov->set_phenotype_vector(phenotypes_);
				  val = call_method(pt_.method, bootstrapped_, *pt_.cov, phenotypes_, pt_.tp, ts);
				  if (val >= original) {
					successes++;
				  }
				}
				// Burn in 10 reps, then start checking the p-value confidence interval
				p = (successes + 2.) / (n + 4.); // Wilson Score estimate
				z = 1.96;
				z2 = z * z;
				sp = (p + z2 / (2 * n)) / (1. + z2 / n);
				ci = z / (1. + z2 / n) * std::sqrt(p * (1. - p) / n + z2 / (4 * n * n)); // Wilson score interval
			  } while (n < 10 || (a <= sp + ci && a >= sp - ci));

			  if (sp + ci < a) {
				success_map_[ts][k]++;
			  }
			}
		  }
		}
	  }
	}
  }
  for (const auto &ts : pt_.gene.get_transcripts()) {
	for (arma::uword i = 0; i < pt_.tp.alpha.n_elem; i++) {
	  for (arma::uword j = 0; j < pt_.tp.ncases.size(); j++) {
		PowerRes pr{pt_.gene.gene_name,
					ts,
					pt_.tp.method,
					pt_.tp.ncases[j],
					pt_.tp.ncontrols[j],
					static_cast<double>(success_map_[ts][j](i)),
					static_cast<double>(pt_.nreps),
					static_cast<double>(success_map_[ts][j](i)) / pt_.nreps,
					pt_.tp.alpha(i)
		};
		pt_.power_res_.emplace_back(pr);
	  }
	}
  }
  std::cerr << "Finished: " << pt_.gene.gene_name << std::endl;
  done_ = true;
}

arma::sp_mat PowerOp::sample(arma::sp_mat &X, arma::uword ncases, arma::uword ncontrols) {
  arma::mat Xmat(X); // Convert to dense matrix for smarter subsetting
  case_idx_ = arma::uvec(ncases, arma::fill::zeros);
  control_idx_ = arma::uvec(ncontrols, arma::fill::zeros);

  // Generates on closed interval [a, b]
  std::uniform_int_distribution<> case_dis(0, cases_.n_elem - 1);
  std::uniform_int_distribution<> control_dis(0, controls_.n_elem - 1);

  for (arma::uword i = 0; i < ncases; i++) {
	case_idx_[i] = static_cast<arma::uword>(case_dis(gen_));
  }
  for (arma::uword i = 0; i < ncontrols; i++) {
	control_idx_[i] = static_cast<arma::uword>(control_dis(gen_));
  }

  arma::uvec uniq_case = arma::unique(case_idx_);
  arma::uvec uniq_control = arma::unique(control_idx_);

  arma::sp_mat ret(arma::join_cols(Xmat.rows(cases_(case_idx_)), Xmat.rows(controls_(control_idx_))));

  return ret;
}

double PowerOp::call_method(Methods &method,
							Gene &gene,
							Covariates &cov,
							arma::vec &phenotypes,
							TaskParams &tp,
							const std::string &k) {
  bool shuffle = false; // Not needed here
  bool detail = false;
  if (tp.method == "BURDEN") {
	return method.BURDEN(gene, k, phenotypes, tp.a, tp.b);
  } else if (tp.method == "CALPHA") {
	return method.CALPHA(gene, phenotypes, k);
  } else if (tp.method == "CMC") {
	return method.CMC(gene, phenotypes, k, tp.maf);
  } else if (tp.method == "RVT1") {
	return method.RVT1(gene, phenotypes, cov.get_covariate_matrix(), cov.get_coef(), k, tp.linear);
  } else if (tp.method == "RVT2") {
	return method.RVT2(gene, phenotypes, cov.get_covariate_matrix(), cov.get_coef(), k, tp.linear);
  } else if (tp.method == "SKAT") {
	return method.SKAT(gene, k, phenotypes, tp.a, tp.b, detail, tp.linear, false, false);
  } else if (tp.method == "SKATO") {
	return method.SKATO(gene, k, phenotypes, tp.a, tp.b, detail, tp.linear);
  } else if (tp.method == "VAAST") {
	return method
		.VAAST(gene,
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
