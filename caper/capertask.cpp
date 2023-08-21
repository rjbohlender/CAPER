//
// Created by Bohlender,Ryan James on 8/22/18.
//

#include "capertask.hpp"
#include "../link/binomial.hpp"
#include "../statistics/bayesianglm.hpp"

CAPERTask::CAPERTask(Stage stage_, Gene gene_,
                     const std::shared_ptr<Covariates> &cov_, TaskParams tp_,
                     arma::uword succ_thresh_, arma::uword nperm_,
                     arma::uword offset_, arma::uword termination_,
                     std::vector<std::vector<int8_t>> &perm_)
    : stage(stage_), gene(std::move(gene_)), cov(cov_), methods(tp_, cov_),
      success_threshold(succ_thresh_), npermutations(nperm_), offset(offset_),
      termination(termination_), permutations(perm_), tp(tp_) {
  for (const auto &ts : gene.get_transcripts()) {
    results[ts] =
        Result(gene.gene_name, ts, !gene.is_polymorphic(ts));
    results[ts].output_stats = tp.output_stats;
    if (tp.seed) {
      permute[ts] = Permute(*tp.seed);
    } else {
      permute[ts] = Permute();
    }
  }
}

CAPERTask::CAPERTask(Stage stage_, Gene &gene_,
                     std::shared_ptr<Covariates> cov_, TaskParams tp_,
                     std::vector<std::vector<int8_t>> &perm_)
    : stage(stage_), gene(std::move(gene_)), cov(cov_), methods(tp_, cov_),
      npermutations(tp_.nperm), permutations(perm_),
      success_threshold(tp_.success_threshold), offset(0),
      termination(tp_.nperm), tp(tp_) {
  for (const auto &transcript : gene.get_transcripts()) {
    results.emplace(
        transcript, Result(gene.gene_name, transcript,
                                          !gene.is_polymorphic(transcript)));
    results[transcript].output_stats = tp.output_stats;
    if (tp.seed) {
      permute[transcript] = Permute(*tp.seed);
    } else {
      permute[transcript] = Permute();
    }
  }
}

Covariates &CAPERTask::get_cov() { return *cov; }

int CAPERTask::max_permutations() {
  auto res = std::max_element(
      results.cbegin(), results.cend(), [](const auto &v1, const auto &v2) {
        return v1.second.permutations < v2.second.permutations;
      });
  return (*res).second.permutations;
}

std::vector<std::vector<int8_t>> &CAPERTask::get_permutations() {
  return permutations;
}

void CAPERTask::cleanup() {
  if (!tp.qtl) {
    Binomial link("logit");
    for (const auto &ts : gene.get_transcripts()) {
      // arma::mat X = arma::mat(arma::sum(gene.genotypes[ts].t(), 0).t());
      // arma::mat D = arma::join_horiz(cov->get_covariate_matrix(), X);
      // BayesianGLM<Binomial> fit(D, cov->get_original_phenotypes(), link);
      // results[ts].odds = std::exp(fit.beta_(fit.beta_.n_elem - 1));
    }
  } else {
    for (const auto &ts : gene.get_transcripts()) {
      // results[ts].odds = 1; // default
    }
  }
  gene.clear(*cov, results, tp);
  // cov_->clear();
  methods.clear(gene.get_transcripts());
}

void CAPERTask::calc_multitranscript_pvalues() {
  // Skip MGIT if only a single transcript
  if (!has_multiple_transcripts()) {
    for (auto &v : results) {
      v.second.mgit_p = v.second.empirical_p;
      v.second.mgit_successes = v.second.successes;
    }
    return;
  }

  // Shorthand
  std::vector<std::string> &transcripts = gene.get_transcripts();

  unsigned long n = transcripts.size();

  int i, j, k;
  double successes = 0;
  double midp_successes = 0;

  arma::mat mgit_pval_mat = arma::mat(max_permutations() + 1, n);

  // For each transcript
  for (i = 0; i < n; i++) {
    // For each statistic
    const std::string &ts = transcripts[i];
    if (!gene.is_polymorphic(ts)) {
      continue;
    }
    int m = static_cast<int>(results[ts].permuted.size()); // Total permutations

    assert(m == max_permutations()); // Sanity check - All equal

    // Append original
    results[ts].permuted.push_back(results[ts].original);
    arma::vec permuted = arma::conv_to<arma::vec>::from(results[ts].permuted);

    // Remove original
    results[ts].permuted.pop_back();

    arma::vec pvals;
    if (tp.analytic) { // Analytic methods return p-values, so ordering needs to
                       // be reversed
      pvals = rank(permuted, "ascend");
    } else {
      pvals = rank(permuted, "descend");
    }

    pvals /= permuted.n_rows;

    try {
      mgit_pval_mat.col(i) = pvals;
    } catch (const std::logic_error &e) {
      std::cerr << "n_row: " << mgit_pval_mat.n_rows
                << " n_col: " << mgit_pval_mat.n_cols << "\n";
      throw(e);
    }
  }

  // Min value of each row
  arma::vec mgit_pval_dist_ = arma::min(mgit_pval_mat, 1);

  for (i = 0; i < n; i++) {
    const std::string &ts = transcripts[i];
    unsigned long m = mgit_pval_dist_.n_rows; // Total permutations

    successes =
        arma::find(mgit_pval_dist_ <= results[ts].empirical_p).eval().n_rows;
    midp_successes =
        arma::find(mgit_pval_dist_ < results[ts].empirical_p).eval().n_rows;
    midp_successes +=
        arma::find(mgit_pval_dist_ == results[ts].empirical_p).eval().n_rows / 2.;

    // Store multi-transcript p-value
    results[ts].mgit_p = (1.0 + successes) / (1.0 + m);
    results[ts].mgit_successes = static_cast<int>(successes);
    results[ts].mgit_midp = (0.5 + midp_successes) / (1.0 + m);
    results[ts].mgit_midp_successes = static_cast<int>(midp_successes);
  }
}

bool CAPERTask::has_multiple_transcripts() {
  return gene.get_transcripts().size() > 1;
}
