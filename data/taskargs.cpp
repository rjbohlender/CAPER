//
// Created by Bohlender,Ryan James on 8/22/18.
//

#include "taskargs.hpp"

TaskArgs::TaskArgs(Stage stage,
				   Gene gene,
				   Covariates cov,
				   TaskParams &tp,
				   int succ_thresh,
				   int s1_perm,
				   int s2_perm,
				   std::vector<std::vector<int32_t>> &perm)
	: stage_(stage),
	  gene_(std::move(gene)),
	  cov_(cov),
	  method_(tp.method, tp.kernel, cov),
	  stage_1_permutations_(s1_perm),
	  stage_2_permutations_(s2_perm),
	  permutations(perm),
	  success_threshold(succ_thresh),
	  stop_check_threshold(succ_thresh),
	  adjust(tp.adjust),
	  tp_(tp) {
  for (const auto &k : gene_.get_transcripts()) {
	results[k] = Result(gene_.get_gene(), k);
	permute[k] = Permute();
  }
}

TaskArgs::TaskArgs(Stage stage,
				   Gene gene,
				   Covariates cov,
				   TaskParams &tp,
				   std::vector<std::vector<int32_t>> &perm)
	: stage_(stage),
	  gene_(std::move(gene)),
	  cov_(cov),
	  method_(tp.method, tp.kernel, cov),
	  stage_1_permutations_(tp.stage_1_permutations),
	  stage_2_permutations_(tp.stage_2_permutations),
	  permutations(perm),
	  success_threshold(tp.success_threshold),
	  stop_check_threshold(tp.success_threshold),
	  adjust(tp.adjust),
	  tp_(tp) {
  for (const auto &k : gene_.get_transcripts()) {
	results[k] = Result(gene_.get_gene(), k);
	permute[k] = Permute();
  }
}
TaskArgs::TaskArgs(const TaskArgs &ta)
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

TaskArgs::TaskArgs(TaskArgs &&ta) noexcept
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

TaskArgs &TaskArgs::operator=(const TaskArgs &rhs) {
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

const Stage &TaskArgs::get_stage() const {
  return stage_;
}

void TaskArgs::set_stage(Stage stage) {
  this->stage_ = stage;
}

Gene &TaskArgs::get_gene() {
  return gene_;
}

Covariates &TaskArgs::get_cov() {
  return cov_;
}

Methods &TaskArgs::get_methods() {
  return method_;
}

int TaskArgs::get_max_permutations() {
  auto res = std::max_element(results.cbegin(),
							  results.cend(),
							  [](const auto &v1, const auto &v2) {
								return v1.second.permutations < v2.second.permutations;
							  });
  return (*res).second.permutations;
}

int TaskArgs::get_npermutations() {
  if (stage_ == Stage::Stage1) {
	return stage_1_permutations_;
  } else if (stage_ == Stage::Stage2) {
	return stage_2_permutations_ - stage_1_permutations_;
  } else {
	return 0;
  }
}

std::vector<std::vector<int32_t>> &TaskArgs::get_permutations() {
  return permutations;
}

void TaskArgs::cleanup() {
  gene_.clear(cov_, results);
  cov_.clear();
  method_.clear(gene_.get_transcripts());
}

Result &TaskArgs::max_original_statistic() {
  auto res = std::max_element(
	  results.begin(),
	  results.end(),
	  [](const auto &res1, const auto &res2) {
		return res1.second.original < res2.second.original;
	  });
  return (*res).second;
}

Result &TaskArgs::min_empirical_pvalue() {
  auto res = std::min_element(
	  results.begin(),
	  results.end(),
	  [](const auto &res1, const auto &res2) {
		return res1.second.empirical_p < res2.second.empirical_p;
	  });
  return (*res).second;
}

Result &TaskArgs::min_mgit_p() {
  auto res = std::min_element(
	  results.begin(),
	  results.end(),
	  [](const auto &res1, const auto &res2) {
		return res1.second.mgit_p < res2.second.mgit_p;
	  });
  return (*res).second;
}

Result &TaskArgs::min_transcript_permutations() {
  auto res = std::min_element(
	  results.begin(),
	  results.end(),
	  [](const auto &res1, const auto &res2) {
		return res1.second.permutations < res2.second.permutations;
	  });
  return (*res).second;
}

Result &TaskArgs::min_transcript_successes() {
  auto res = std::min_element(
	  results.begin(),
	  results.end(),
	  [](const auto &res1, const auto &res2) {
		return res1.second.successes < res2.second.successes;
	  });
  return (*res).second;
}

double TaskArgs::get_mgit_pvalue(const std::string &k) {
  if (has_multiple_transcripts()) {
	return results[k].mgit_p;
  } else {
	return results[k].empirical_p;
  }
}

int TaskArgs::get_remaining() {
  if (stage_ == Stage::Stage1) {
	return stage_2_permutations_ > 0 ? stage_2_permutations_ - stage_1_permutations_ : 0;
  } else {
	return 0;
  }
}

void TaskArgs::calc_multitranscript_pvalues() {
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
	int m = static_cast<int>(results[ts].permuted.size());  // Total permutations

	assert(m == get_max_permutations()); // Sanity check - All equal

	// Append original
	results[ts].permuted.push_back(results[ts].original);
	arma::vec permuted = arma::conv_to<arma::vec>::from(results[ts].permuted);

	// Remove original
	results[ts].permuted.pop_back();

	arma::vec pvals;
	if (method_.str() == "SKATO" || method_.str() == "SKAT") {
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

bool TaskArgs::has_multiple_transcripts() {
  return gene_.get_transcripts().size() > 1;
}

Permute &TaskArgs::get_permute(const std::string &k) {
  return permute[k];
}

int TaskArgs::get_a() {
  return tp_.a;
}

int TaskArgs::get_b() {
  return tp_.b;
}


TaskParams &TaskArgs::get_tp() {
  return tp_;
}
