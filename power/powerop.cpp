//
// Created by Bohlender,Ryan James on 2019-03-19.
//

#include "powerop.hpp"
#include <iostream>
#include <iomanip>

PowerOp::PowerOp(PowerTask &pt, std::shared_ptr<PowerReporter> reporter, double seed, bool verbose)
	: caperTask(pt), reporter_(reporter), gen_(seed), bootstrapped_(caperTask.gene), done_(false) {
}

auto PowerOp::run() -> void {
  power();
}

auto PowerOp::finish() -> void {
  // cleanup
  if (!caperTask.tp.gene_list) {
    reporter_->sync_write_power(caperTask.power_res_);
  }
}

auto PowerOp::is_done() const -> bool {
  return done_;
}

auto PowerOp::get_task() -> PowerTask {
  return caperTask;
}

auto PowerOp::power() -> void {
  if (cases_.empty()) {
	cases_ = arma::find(caperTask.cov->get_original_phenotypes() > 0);
  }
  if (controls_.empty()) {
	controls_ = arma::find(caperTask.cov->get_original_phenotypes() == 0);
  }

  for (const auto &ts : caperTask.gene.get_transcripts()) {
        const bool method_returns_pvalue =
            caperTask.tp.method == "SKATO" || caperTask.tp.method == "SKAT" ||
            caperTask.tp.method == "CMC" || caperTask.tp.method == "RVT1" ||
            caperTask.tp.method == "RVT2";
        if (caperTask.tp.alpha.size() > 1) {
	  for (arma::uword k = 0; k < caperTask.tp.ncases.size(); k++) {
		// Alpha vector to track the successes
		success_map_[ts].push_back(arma::vec(caperTask.tp.alpha.size(), arma::fill::zeros));
		arma::uword ncases = caperTask.tp.ncases[k];
		arma::uword ncontrols = caperTask.tp.ncontrols[k];

		phenotypes_.zeros(ncases + ncontrols);
		phenotypes_(arma::span(0, ncases - 1)).fill(1);

		arma::vec odds = caperTask.cov->get_odds();
		for (arma::uword i = 0; i < caperTask.nreps; i++) {
		  bootstrapped_.set_matrix(ts, sample(caperTask.gene.genotypes[ts], ncases, ncontrols));
		  // Original for this bootstrap replicate
		  double original = call_method(caperTask.method, bootstrapped_, *caperTask.cov, phenotypes_,
                                                caperTask.tp, ts);
		  // Need to get a p-value via permutation for some methods
		  // If the method returns a p-value
                  if (method_returns_pvalue) {
                        arma::uvec passing = arma::find(caperTask.tp.alpha >= original);
                        if (!passing.is_empty()) {
                          success_map_[ts][k](passing) += 1;
                        }
                  } else {
			// We need to permute to get a p-value
			// Cease permutation if the p-value ci excludes alpha above or below
			// If below, call it a success
			int block = static_cast<int>(
                            caperTask.tp.nperm); // Permutation block size
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
                                caperTask.tp.bin_epsilon);

			for (arma::uword j = 0; j < block; j++) {
			  phenotypes_ = arma::conv_to<arma::vec>::from(permutations[j]);
                          caperTask.cov->set_phenotype_vector(phenotypes_);
			  val = call_method(caperTask.method, bootstrapped_, *caperTask.cov, phenotypes_,
                                            caperTask.tp, ts);
			  if (val >= original) {
				successes++;
			  }
			}

			p = (successes + 2.) / (block + 4.); // Wilson Score estimate
			success_map_[ts][k](arma::find(caperTask.tp.alpha >= p)) += 1;
		  }
		}
	  }
        } else {
          if (caperTask.tp.alpha.n_elem == 0) {
                continue;
          }
          double alpha_threshold = caperTask.tp.alpha(0);
          for (arma::uword k = 0; k < caperTask.tp.ncases.size(); k++) {
                success_map_[ts].push_back(
                    arma::vec(caperTask.tp.alpha.n_elem, arma::fill::zeros));
                arma::uword ncases = caperTask.tp.ncases[k];
                arma::uword ncontrols = caperTask.tp.ncontrols[k];

                phenotypes_.zeros(ncases + ncontrols);
                phenotypes_(arma::span(0, ncases - 1)).fill(1);

                arma::vec odds = caperTask.cov->get_odds();
                for (arma::uword i = 0; i < caperTask.nreps; i++) {
                  bootstrapped_.set_matrix(
                      ts, sample(caperTask.gene.genotypes[ts], ncases, ncontrols));
                  double original =
                      call_method(caperTask.method, bootstrapped_, *caperTask.cov,
                                  phenotypes_, caperTask.tp, ts);
                  if (method_returns_pvalue) {
                        arma::uvec passing =
                            arma::find(caperTask.tp.alpha >= original);
                        if (!passing.is_empty()) {
                          success_map_[ts][k](passing) += 1;
                        }
                  } else {
                        arma::uword n = 0;
                        int block = 5; // Permutation block size
                        double successes = 0;
                        double p, sp;
                        double z, z2;
                        double val;
                        double ci;

                        arma::uvec idx =
                            arma::join_vert(cases_(case_idx_), controls_(control_idx_));
                        arma::vec bootstrapped_odds = odds(idx);
                        do {
                          std::vector<std::vector<int8_t>> permutations =
                              permute_.epsilon_permutation(
                                  block, bootstrapped_odds, ncases, ts,
                                  caperTask.tp.bin_epsilon);

                          for (arma::uword j = 0; j < block; j++) {
                                n++;
                                phenotypes_ =
                                    arma::conv_to<arma::vec>::from(permutations[j]);
                                caperTask.cov->set_phenotype_vector(phenotypes_);
                                val = call_method(caperTask.method, bootstrapped_,
                                                  *caperTask.cov, phenotypes_,
                                                  caperTask.tp, ts);
                                if (val >= original) {
                                  successes++;
                                }
                          }
                          p = (successes + 2.) / (n + 4.); // Wilson Score estimate
                          z = 1.96;
                          z2 = z * z;
                          sp = (p + z2 / (2 * n)) / (1. + z2 / n);
                          ci = z / (1. + z2 / n) *
                               std::sqrt(p * (1. - p) / n + z2 / (4 * n * n));
                        } while (n < 10 || (alpha_threshold <= sp + ci &&
                                            alpha_threshold >= sp - ci));

                        if (sp + ci < alpha_threshold) {
                          success_map_[ts][k](0) += 1;
                        }
                  }
                }
          }
        }
  }
  for (const auto &ts : caperTask.gene.get_transcripts()) {
	for (arma::uword i = 0; i < caperTask.tp.alpha.n_elem; i++) {
	  for (arma::uword j = 0; j < caperTask.tp.ncases.size(); j++) {
		PowerRes pr{caperTask.gene.gene_name,
					ts,
                        caperTask.tp.method,
                        caperTask.tp.ncases[j],
                        caperTask.tp.ncontrols[j],
					static_cast<double>(success_map_[ts][j](i)),
					static_cast<double>(caperTask.nreps),
					static_cast<double>(success_map_[ts][j](i)) /
                            caperTask.nreps,
                        caperTask.tp.alpha(i)
		};
                caperTask.power_res_.emplace_back(pr);
	  }
	}
  }
  std::cerr << "Finished: " << caperTask.gene.gene_name << std::endl;
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
	return method.BURDEN(gene, phenotypes, k, false);
  } else if (tp.method == "CALPHA") {
	return method.CALPHA(gene, phenotypes, k);
  } else if (tp.method == "CMC") {
	return method.CMC(gene, phenotypes, k, tp.maf);
  } else if (tp.method == "RVT1") {
	return method.RVT1(gene, phenotypes, cov.get_covariate_matrix(), cov.get_coef(), k, tp.qtl);
  } else if (tp.method == "RVT2") {
	return method.RVT2(gene, phenotypes, cov.get_covariate_matrix(), cov.get_coef(), k, tp.qtl);
  } else if (tp.method == "SKAT") {
	return method.SKAT(gene, phenotypes, k, tp.a, tp.b, detail, tp.qtl,
                           false);
  } else if (tp.method == "SKATO") {
	return method.SKATO(gene, phenotypes, k, tp.a, tp.b, detail, tp.qtl);
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
			   tp.alternate_grouping);
  } else if (tp.method == "VT") {
	return method.VT(gene, phenotypes, k);
  } else if (tp.method == "WSS") {
	return method.WSS(gene, phenotypes, k);
  } else {
	throw (std::logic_error("Failed to find method."));
  }
}
