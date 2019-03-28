//
// Created by Bohlender,Ryan James on 2019-03-19.
//

#include "power.hpp"
#include <iostream>
#include <iomanip>

Power::Power(Methods &method, Gene &gene, Covariates &cov, TaskParams &tp, arma::uword nreps)
: tp_(tp), gen_(std::random_device{}()), bootstrapped_(gene), cov_(cov) {
  if(cases_.empty()) {
    cases_ = arma::find(cov_.get_original_phenotypes() > 0);
  }
  if(controls_.empty()) {
    controls_ = arma::find(cov_.get_original_phenotypes() == 0);
  }

  for(const auto &ts : gene.get_transcripts()) {
    for(const auto &a : tp.alpha) {
      for(arma::uword k = 0; k < tp.ncases.size(); k++) {
        arma::uword ncases = tp.ncases[k];
        arma::uword ncontrols = tp.ncontrols[k];

        phenotypes_.zeros(ncases + ncontrols);
        phenotypes_(arma::span(0, ncases - 1)).fill(1);

        arma::vec odds = cov_.get_odds();
        // Reset map for this set of bootstraps
        success_map_[ts] = 0;
        for(arma::uword i = 0; i < nreps; i++) {
          bootstrapped_.set_matrix(ts, sample(gene.get_matrix(ts), ncases, ncontrols));
          // Original for this bootstrap replicate
          double original = call_method(method, bootstrapped_, cov_, phenotypes_, tp, ts);
          // Need to get a p-value via permutation for some methods
          // If the method returns a p-value
          if (tp_.method == "SKATO" || tp_.method == "SKAT" || tp_.method == "CMC" || tp_.method == "RVT1" || tp_.method == "RVT2") {
            double val = call_method(method, bootstrapped_, cov_, phenotypes_, tp, ts);
            if(val < original) {
              success_map_[ts]++;
            }
          } else {
            // We need to permute to get a p-value
            // Cease permutation if the p-value ci excludes alpha above or below
            // If below, call it a success
            arma::uword n = 0;
            int block = 5; // Permutation block size
            double successes = 0;
            double p;
            double val;
            double ci;

            arma::uvec idx = arma::join_vert(cases_(case_idx_), controls_(control_idx_));
            arma::vec bootstrapped_odds = odds(idx);
            // Indices of minor allele carriers for grouping on non-carriers
            arma::uvec mac_idx = arma::find(arma::sum(arma::mat(bootstrapped_.get_matrix(ts)), 1) > 0);
            arma::uvec maj_idx = arma::find(arma::sum(arma::mat(bootstrapped_.get_matrix(ts)), 1) == 0);
            do {
              // Shuffle phenotypes
              std::vector<std::vector<int32_t>> permutations = permute_.permutations_maj_bin(
                  block,
                  bootstrapped_odds,
                  ncases,
                  mac_idx,
                  maj_idx);

              for (arma::uword j = 0; j < block; j++) {
                n++;
                phenotypes_ = arma::conv_to<arma::vec>::from(permutations[j]);
                cov_.set_phenotype_vector(phenotypes_);
                val = call_method(method, bootstrapped_, cov_, phenotypes_, tp, ts);
                if (val >= original) {
                  successes++;
                }
              }
              // Burn in 10 reps, then start checking the p-value confidence interval
              p = (successes + 1.) / (n + 1.);
              ci = 1.96 * std::sqrt(p * (1. - p) / n);
            } while(n < 10 || (a <= p + ci && a >= p - ci));

            if (p + ci < a) {
              success_map_[ts]++;
            }
          }
        }
        for (const auto &k : gene.get_transcripts()) {
          PowerRes pr {
              gene.get_gene(),
              k,
              tp_.method,
              ncases,
              ncontrols,
              static_cast<double>(success_map_[k]),
              static_cast<double>(nreps),
              static_cast<double>(success_map_[k]) / nreps,
              a
          };
          power_res_.emplace_back(pr);
        }
      }
    }
  }
  std::cerr << "Finished: " << gene.get_gene() << std::endl;
}

arma::sp_mat Power::sample(arma::sp_mat &X, arma::uword ncases, arma::uword ncontrols) {
  arma::mat Xmat(X); // Convert to dense matrix for smarter subsetting
  case_idx_ = arma::uvec(ncases, arma::fill::zeros);
  control_idx_ = arma::uvec(ncontrols, arma::fill::zeros);

  // Generates on closed interval [a, b]
  std::uniform_int_distribution<> case_dis(0, cases_.n_elem - 1);
  std::uniform_int_distribution<> control_dis(0, controls_.n_elem - 1);

  for(arma::uword i = 0; i < ncases; i++) {
    case_idx_[i] = static_cast<arma::uword>(case_dis(gen_));
  }
  for(arma::uword i = 0; i < ncontrols; i++) {
    control_idx_[i] = static_cast<arma::uword>(control_dis(gen_));
  }

  arma::uvec uniq_case = arma::unique(case_idx_);
  arma::uvec uniq_control = arma::unique(control_idx_);

  arma::sp_mat ret(arma::join_cols(Xmat.rows(cases_(case_idx_)), Xmat.rows(controls_(control_idx_))));

  return ret;
}
double Power::call_method(Methods &method,
                          Gene &gene,
                          Covariates &cov,
                          arma::vec &phenotypes,
                          TaskParams &tp,
                          const std::string &k) {
  bool shuffle = false; // Not needed here
  bool detail = false;
  if (tp.method == "BURDEN") {
    return method.BURDEN(gene, k, shuffle, tp.a, tp.b);
  } else if (tp.method == "CALPHA") {
    return method.CALPHA(gene, phenotypes, k);
  } else if (tp.method == "CMC") {
    return method.CMC(gene, phenotypes, k, tp.maf);
  } else if (tp.method == "RVT1") {
    return method.RVT1(gene, phenotypes, cov.get_covariate_matrix(), k, tp.linear);
  } else if (tp.method == "RVT2") {
    return method.RVT2(gene, phenotypes, cov.get_covariate_matrix(), k, tp.linear);
  } else if (tp.method == "SKAT") {
    return method.SKATR(gene, k, shuffle, tp.a, tp.b, detail, tp.linear, false);
  } else if (tp.method == "SKATO") {
    return method.SKATRO(gene, k, shuffle, tp.a, tp.b, detail, tp.linear);
  } else if (tp.method == "VAAST") {
    return method.Vaast(gene, phenotypes, k, tp.score_only_minor, tp.score_only_alternative, 2.0, tp.group_size, detail);
  } else if (tp.method == "VT") {
    return method.VT(gene, k, shuffle);
  } else if (tp.method == "WSS") {
    return method.WSS(gene, phenotypes, k);
  } else {
    throw(std::logic_error("Failed to find method."));
  }
}

std::vector<PowerRes> Power::get_results() {
  return power_res_;
}
