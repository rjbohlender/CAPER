//
// Created by Bohlender,Ryan James on 8/22/18.
//

#include "carvatask.hpp"
#include "../statistics/bayesianglm.hpp"
#include "../link/binomial.hpp"

const std::set<std::string> CARVATask::pvalue_methods_{
	"SKAT",
	"SKATO",
	"CMC"
};

CARVATask::CARVATask(Stage stage,
					 Gene gene,
					 const std::shared_ptr<Covariates> &cov,
					 TaskParams &tp,
					 arma::uword succ_thresh,
					 arma::uword s1_perm,
					 arma::uword s2_perm,
					 std::vector<std::vector<int32_t>> &perm)
	: stage_(stage),
	  gene_(std::move(gene)),
	  cov_(cov),
	  method_(tp, cov),
	  stage_1_permutations_(s1_perm),
	  stage_2_permutations_(s2_perm),
	  permutations(perm),
	  success_threshold(succ_thresh),
	  stop_check_threshold(succ_thresh),
	  adjust(tp.adjust),
	  tp_(tp) {
  for (const auto &k : gene_.get_transcripts()) {
	results[k] = Result(gene_.get_gene(), k, !gene_.is_polymorphic(k));
	results[k].output_stats = tp.output_stats;
	permute[k] = Permute();
  }
}

CARVATask::CARVATask(Stage stage,
					 Gene &gene,
					 const std::shared_ptr<Covariates> &cov,
					 TaskParams &tp,
					 std::vector<std::vector<int32_t>> &perm)
	: stage_(stage),
	  gene_(std::move(gene)),
	  cov_(cov),
	  method_(tp, cov),
	  stage_1_permutations_(tp.stage_1_permutations),
	  stage_2_permutations_(tp.stage_2_permutations),
	  permutations(perm),
	  success_threshold(tp.success_threshold),
	  stop_check_threshold(tp.success_threshold),
	  adjust(tp.adjust),
	  tp_(tp) {
  for (const auto &k : gene_.get_transcripts()) {
	results.emplace(std::make_pair(k, Result(gene_.get_gene(), k, !gene_.is_polymorphic(k))));
	results[k].output_stats = tp.output_stats;
	permute[k] = Permute();
  }
}
CARVATask::CARVATask(const CARVATask &ta)
	: results(ta.results),
	  permute(ta.permute),
	  stage_(ta.stage_),
	  gene_(ta.gene_),
	  cov_(ta.cov_),
	  method_(ta.method_),
	  stage_1_permutations_(ta.stage_1_permutations_),
	  stage_2_permutations_(ta.stage_2_permutations_),
	  permutations(ta.permutations),
	  success_threshold(ta.success_threshold),
	  stop_check_threshold(ta.stop_check_threshold),
	  adjust(ta.adjust),
	  tp_(ta.tp_) {}

CARVATask::CARVATask(CARVATask &&ta) noexcept
	: results(std::move(ta.results)),
	  permute(std::move(ta.permute)),
	  stage_(ta.stage_),
	  gene_(std::move(ta.gene_)),
	  cov_(ta.cov_),
	  method_(std::move(ta.method_)),
	  stage_1_permutations_(ta.stage_1_permutations_),
	  stage_2_permutations_(ta.stage_2_permutations_),
	  permutations(ta.permutations),
	  success_threshold(ta.success_threshold),
	  stop_check_threshold(ta.stop_check_threshold),
	  adjust(ta.adjust),
	  tp_(ta.tp_) {}

CARVATask &CARVATask::operator=(const CARVATask &rhs) {
  stage_ = rhs.stage_;
  gene_ = rhs.gene_;
  cov_ = rhs.cov_;
  results = rhs.results;
  permute = rhs.permute;
  stage_1_permutations_ = rhs.stage_1_permutations_;
  stage_2_permutations_ = rhs.stage_2_permutations_;
  permutations = rhs.permutations;
  success_threshold = rhs.success_threshold;
  stop_check_threshold = rhs.stop_check_threshold;
  adjust = rhs.adjust;
  tp_ = rhs.tp_;

  return *this;
}

const Stage &CARVATask::get_stage() const {
  return stage_;
}

void CARVATask::set_stage(Stage stage) {
  this->stage_ = stage;
}

Gene &CARVATask::get_gene() {
  return gene_;
}

Covariates &CARVATask::get_cov() {
  return *cov_;
}

Methods &CARVATask::get_methods() {
  return method_;
}

int CARVATask::get_max_permutations() {
  auto res = std::max_element(results.cbegin(),
							  results.cend(),
							  [](const auto &v1, const auto &v2) {
								return v1.second.permutations < v2.second.permutations;
							  });
  return (*res).second.permutations;
}

int CARVATask::get_npermutations() {
  if (stage_ == Stage::Stage1) {
	return stage_1_permutations_;
  } else if (stage_ == Stage::Stage2) {
	return stage_2_permutations_ - stage_1_permutations_;
  } else {
	return 0;
  }
}

std::vector<std::vector<int32_t>> &CARVATask::get_permutations() {
  return permutations;
}

void CARVATask::cleanup() {
  if (!tp_.linear) {
	Binomial link("logit");
	for (const auto &ts : gene_.get_transcripts()) {
	  arma::mat X = arma::mat(arma::sum(gene_.get_matrix(ts).t(), 0).t());
	  arma::mat D = arma::join_horiz(cov_->get_covariate_matrix(), X);
	  BayesianGLM<Binomial> fit(D, cov_->get_original_phenotypes(), link);
	  results[ts].set_odds(std::exp(fit.beta_(fit.beta_.n_elem - 1)));
	}
  } else {
	for (const auto &ts : gene_.get_transcripts()) {
	  results[ts].set_odds(1); // default
	}
  }
  gene_.clear(*cov_, results, tp_);
  //cov_->clear();
  method_.clear(gene_.get_transcripts());
}

Result &CARVATask::max_original_statistic() {
  auto res = std::max_element(
	  results.begin(),
	  results.end(),
	  [](const auto &res1, const auto &res2) {
		return res1.second.original < res2.second.original;
	  });
  return (*res).second;
}

Result &CARVATask::min_empirical_pvalue() {
  auto res = std::min_element(
	  results.begin(),
	  results.end(),
	  [](const auto &res1, const auto &res2) {
		return res1.second.empirical_p < res2.second.empirical_p;
	  });
  return (*res).second;
}

Result &CARVATask::min_mgit_p() {
  auto res = std::min_element(
	  results.begin(),
	  results.end(),
	  [](const auto &res1, const auto &res2) {
		return res1.second.mgit_p < res2.second.mgit_p;
	  });
  return (*res).second;
}

Result &CARVATask::min_transcript_permutations() {
  auto res = std::min_element(
	  results.begin(),
	  results.end(),
	  [](const auto &res1, const auto &res2) {
		return res1.second.permutations < res2.second.permutations;
	  });
  return (*res).second;
}

Result &CARVATask::min_transcript_successes() {
  auto res = std::min_element(
	  results.begin(),
	  results.end(),
	  [](const auto &res1, const auto &res2) {
		return res1.second.successes < res2.second.successes;
	  });
  return (*res).second;
}

double CARVATask::get_mgit_pvalue(const std::string &k) {
  if (has_multiple_transcripts()) {
	return results[k].mgit_p;
  } else {
	return results[k].empirical_p;
  }
}

int CARVATask::get_remaining() {
  if (stage_ == Stage::Stage1) {
	return stage_2_permutations_ > 0 ? stage_2_permutations_ - stage_1_permutations_ : 0;
  } else {
	return 0;
  }
}

void CARVATask::calc_multitranscript_pvalues() {
  // Skip MGIT if only a single transcript
  if (!has_multiple_transcripts()) {
	for (auto &v : results) {
	  v.second.mgit_p = v.second.empirical_p;
	  v.second.mgit_successes = v.second.successes;
	}
	return;
  }

  // Shorthand
  std::vector<std::string> &transcripts = gene_.get_transcripts();

  unsigned long n = transcripts.size();

  int i, j, k;
  double successes;

  arma::mat mgit_pval_mat = arma::mat(get_max_permutations() + 1, n);

  // For each transcript
  for (i = 0; i < n; i++) {
	// For each statistic
	const std::string &ts = transcripts[i];
	if (!gene_.is_polymorphic(ts)) {
	  continue;
	}
	int m = static_cast<int>(results[ts].permuted.size());  // Total permutations

	assert(m == get_max_permutations()); // Sanity check - All equal

	// Append original
	results[ts].permuted.push_back(results[ts].original);
	arma::vec permuted = arma::conv_to<arma::vec>::from(results[ts].permuted);

	// Remove original
	results[ts].permuted.pop_back();

	arma::vec pvals;
	if (pvalue_methods_.find(method_.str()) != pvalue_methods_.end()) {
	  // SKATO Returns pvalues so reverse success criteria
	  pvals = rank(permuted, "ascend");
	} else {
	  pvals = rank(permuted, "descend");
	}

	pvals /= permuted.n_rows;

	try {
	  mgit_pval_mat.col(i) = pvals;
	} catch (const std::logic_error &e) {
	  std::cerr << "n_row: " << mgit_pval_mat.n_rows << " n_col: " << mgit_pval_mat.n_cols << "\n";
	  throw (e);
	}
  }

  arma::vec mgit_pval_dist_ = arma::min(mgit_pval_mat, 1);

  for (i = 0; i < n; i++) {
	const std::string &ts = transcripts[i];
	unsigned long m = mgit_pval_dist_.n_rows;  // Total permutations

	successes = arma::find(mgit_pval_dist_ <= results[ts].empirical_p).eval().n_rows;

	// Store multi-transcript p-value
	results[ts].mgit_p = (1.0 + successes) / (1.0 + m);
	results[ts].mgit_successes = static_cast<int>(successes);
  }
}

bool CARVATask::has_multiple_transcripts() {
  return gene_.get_transcripts().size() > 1;
}

Permute &CARVATask::get_permute(const std::string &k) {
  return permute[k];
}

TaskParams &CARVATask::get_tp() {
  return tp_;
}
