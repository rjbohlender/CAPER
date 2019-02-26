//
// Created by Bohlender,Ryan James on 9/17/18.
//

#include "main_support.hpp"

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <chrono>
#include <iomanip>
#include <iostream>

#include "../data/bed.hpp"
#include "../data/weight.hpp"

using namespace std::chrono_literals;

void calc_mgit_pvalues(std::map<std::string, Result> &results,
					   std::vector<std::string> &transcripts,
					   const std::string &method) {
  // Skip MGIT if only a single transcript
  if (transcripts.size() == 1) {
	for (const auto &tr : transcripts) {
	  auto &v = results[tr];
	  v.mgit_p = v.empirical_p;
	  v.mgit_successes = v.successes;
	}
	return;
  }

  // Shorthand
  unsigned long n = transcripts.size();
  int max_perm = 0;
  int i, j, k;
  double successes;

  // Get max_perm
  for (const auto &tr : transcripts) {
	if (results[tr].permutations > max_perm) {
	  max_perm = results[tr].permutations;
	}
  }

  arma::mat mgit_pval_mat = arma::mat(max_perm + 1, n);

  // For each transcript
  for (i = 0; i < n; i++) {
	// For each statistic
	const std::string &ts = transcripts[i];
	int m = static_cast<int>(results[ts].permuted.size());  // Total permutations

	assert(m == max_perm); // Sanity check - All equal

	results[ts].permuted.push_back(results[ts].original);
	arma::vec permuted = arma::conv_to<arma::vec>::from(results[ts].permuted);

	arma::vec pvals;
	if (method == "SKATO" || method == "SKAT") {
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

void write_simple(TaskQueue &tq, arma::uword ntranscripts, arma::uword ngenes, TaskParams &tp) {
  // Build output_path
  std::stringstream simple_path_ss;
  simple_path_ss << tp.output_path << "/" << tp.method << ".simple";
  std::ofstream simple_ofs(simple_path_ss.str());

  // Holds unfinished genes
  std::vector<std::string> unfinished;

  arma::uword ncases = tq.get_results()[0].get_cov().get_ncases();
  arma::uword nconts = tq.get_results()[0].get_cov().get_nsamples() - ncases;

  simple_ofs << "##\ttotal_control\t" << nconts << std::endl;
  simple_ofs << "##\ttotal_seq_aff\t" << ncases << std::endl;
  simple_ofs << std::setw(20) << std::left << "Rank";
  simple_ofs << std::setw(20) << "Gene";
  simple_ofs << std::setw(20) << "Transcript";
  simple_ofs << std::setw(20) << "Original";
  simple_ofs << std::setw(20) << "Empirical_P";
  simple_ofs << std::setw(20) << "Empirical_MidP";
  simple_ofs << std::setw(20) << "Successes";
  simple_ofs << std::setw(20) << "Permutations";
  simple_ofs << std::setw(20) << "MGIT";
  simple_ofs << std::setw(20) << "MGIT_Successes" << std::endl;
  if (!tp.gene_list) {
	// Print header and formatted results
	double permutation_mean = 0;
	double permutation_variance = 0;

	std::vector<Result> all_results;
	for (auto &v : tq.get_results()) {
	  // Result &res = v.min_empirical_pvalue();
	  bool unfin = false;
	  for (auto &u : v.results) {
		all_results.emplace_back(std::move(u.second));

		if (u.second.successes < tp.success_threshold && !unfin) {
		  unfinished.push_back(all_results.back().gene);
		  unfin = true;
		}
	  }

	  for (const auto &k : v.get_gene().get_transcripts()) {
		permutation_mean += v.results[k].permutations;
	  }
	}
	permutation_mean /= ntranscripts;
	for (auto &v : tq.get_results()) {
	  for (const auto &k : v.get_gene().get_transcripts()) {
		permutation_variance += std::pow(v.results[k].permutations - permutation_mean, 2) / ntranscripts;
	  }
	}
	// Sort by midp and score, then output
	std::sort(all_results.begin(),
			  all_results.end(),
			  [](Result &a, Result &b) {
	  if(a.empirical_midp != b.empirical_midp) {
		return a.empirical_midp < b.empirical_midp;
	  } else {
	    return b.original < a.original;
	  }
	});
	// Report results
	int rank = 1;
	for(auto &v : all_results) {
	  v.set_rank(rank);
	  if(tp.testable) {
	    if(v.testable) {
		  simple_ofs << v;
	    } else {
		  continue;
	    }
	  } else {
		simple_ofs << v;
	  }
	  rank++;
	}


	std::cerr << "Permutation mean: " << permutation_mean << std::endl;
	std::cerr << "Permutation sd: " << std::sqrt(permutation_variance) << std::endl;
	std::cerr << "Genes submitted: " << ngenes << std::endl;
	std::cerr << "Transcripts submitted: " << ntranscripts << std::endl;
  } else {
	std::map<std::string, Result> res;
	// Merge gene results
	do {
	  auto &v = tq.get_results().front();

	  std::vector<std::string> transcripts = v.get_gene().get_transcripts();

	  for (const auto &k : transcripts) {
		res[k] = std::move(v.results[k]);
	  }

	  // Check for additional results
	  if (tq.get_results().size() > 1) {
		std::vector<int> to_pop;
		for (int i = 1; i < tq.get_results().size(); i++) {
		  auto &cur_transcripts = tq.get_results()[i].get_gene().get_transcripts();
		  if (std::equal(cur_transcripts.cbegin(), cur_transcripts.cend(), transcripts.cbegin())) {
			// Mark for removal
			to_pop.push_back(i);
			// Combine
			for (const auto &k : cur_transcripts) {
			  res[k].combine(tq.get_results()[i].results[k]);
			}
		  }
		}
		// Remove remaining
		for (auto rit = to_pop.rbegin(); rit != to_pop.rend(); rit++) {
		  tq.get_results().erase(tq.get_results().begin() + *rit);
		}
		calc_mgit_pvalues(res, transcripts, v.get_methods().str());
	  }
	  tq.get_results().erase(tq.get_results().begin());

	  for (const auto &k : transcripts) {
		if (res[k].successes < tp.success_threshold) {
		  unfinished.push_back(res[k].gene);
		  break;
		}
	  }
	  // Renew iterators
	} while (!tq.get_results().empty());

	// Get results
	std::vector<Result> all_results;
	all_results.reserve(res.size());
	for (auto &v : res) {
	  all_results.emplace_back(std::move(v.second));
	}
	// Sort results by midp
	std::sort(all_results.begin(),
			  all_results.end(),
			  [](Result &a, Result &b) { return a.empirical_midp < b.empirical_midp; });
	// Output
	int rank = 1;
	for (auto &v : all_results) {
	  v.set_rank(rank);
	  simple_ofs << v;
	  rank++;
	}
  }
  simple_ofs.flush();
  simple_ofs.close();

  if (!unfinished.empty()) {
	// Print command to run unfinished
	std::stringstream uf_ss;
	uf_ss << tp.program_path << " ";
	uf_ss << "-g " << tp.genotypes_path << " ";
	uf_ss << "-c " << tp.covariates_path << " ";
	if (tp.bed) {
	  uf_ss << "-b " << *tp.bed << " ";
	}
	if (tp.weight) {
	  uf_ss << "-w " << *tp.weight << " ";
	}
	uf_ss << "-1 0 "; // Skip stage 1
	uf_ss << "-2 " << tp.total_permutations * 10 << " ";
	uf_ss << "-m " << tp.method << " ";
	uf_ss << "-t " << tp.nthreads << " ";
	if (tp.adjust) {
	  uf_ss << "-n ";
	}
	if (!tp.verbose) {
	  uf_ss << "-q ";
	}
	if (tp.a != 1 || tp.b != 25) {
	  uf_ss << "--beta_weights " << tp.a << "," << tp.b << " ";
	}

	uf_ss << "-l ";
	for (int i = 0; i < unfinished.size(); i++) {
	  if (i == unfinished.size() - 1) {
		uf_ss << unfinished[i];
	  } else {
		uf_ss << unfinished[i] << ",";
	  }
	}
	std::cerr << "Some genes did not reach the success threshold. Run the following command to check those genes."
			  << std::endl;
	std::cerr << uf_ss.str() << std::endl;
  }
}

void write_detail(TaskQueue &tq, TaskParams &tp) {
  std::stringstream detail_path_ss;
  detail_path_ss << tp.output_path << "/" << tp.method << ".detail";

  std::ofstream detail(detail_path_ss.str());

  std::string header =
	  "#Gene\tTranscripts\tVariant\tScore\tAF\tcase_ref\tcase_alt\tcontrol_ref\tcontrol_alt\tcase_list\tcontrol_list";
  detail << header << std::endl;

  if (!tp.gene_list) {
	int i = 0; // For each gene
	for (auto &v : tq.get_results()) {
	  detail << v.get_gene().get_detail();
	  // Print sample / index map at the end
	  if (i == tq.get_results().size() - 1) {
		detail << "## Sample Index Map" << std::endl;
		int j = 0;
		for (auto &s : v.get_gene().get_samples()) {
		  detail << "#\t" << s << "\t" << j << std::endl;
		  j++;
		}
	  }
	  i++;
	}
  } else {
	// Merge gene results
	do {
	  auto &v = tq.get_results_duplicate().front();

	  detail << v.get_gene().get_detail();
	  std::vector<std::string> transcripts = v.get_gene().get_transcripts();

	  // Check for additional results
	  if (tq.get_results_duplicate().size() > 1) {
		std::vector<int> to_pop;
		for (int i = 1; i < tq.get_results_duplicate().size(); i++) {
		  auto &cur_transcripts = tq.get_results_duplicate()[i].get_gene().get_transcripts();
		  if (std::equal(cur_transcripts.cbegin(), cur_transcripts.cend(), transcripts.cbegin())) {
			// Mark for removal
			to_pop.push_back(i);
		  }
		}
		// Remove remaining
		for (auto rit = to_pop.rbegin(); rit != to_pop.rend(); rit++) {
		  tq.get_results_duplicate().erase(tq.get_results_duplicate().begin() + *rit);
		}
	  }
	  // Print before erasing
	  if (tq.get_results_duplicate().size() == 1) {
		detail << "## Sample Index Map" << std::endl;
		int j = 0;
		for (auto &s : v.get_gene().get_samples()) {
		  detail << "#\t" << s << "\t" << j << std::endl;
		  j++;
		}
	  }
	  tq.get_results_duplicate().erase(tq.get_results_duplicate().begin());

	} while (!tq.get_results_duplicate().empty());
  }
}

