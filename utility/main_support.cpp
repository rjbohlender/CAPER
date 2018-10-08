//
// Created by Bohlender,Ryan James on 9/17/18.
//

#include "main_support.hpp"

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
	// Sort by midp and output
	std::sort(all_results.begin(),
			  all_results.end(),
			  [](Result &a, Result &b) { return a.empirical_midp < b.empirical_midp; });
	// Report results
	int rank = 1;
	for(auto &v : all_results) {
	  v.set_rank(rank);
	  simple_ofs << v;
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
	for (const auto &v : res) {
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
	  uf_ss << "-q";
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

void initialize_jobs(TaskParams &tp,
					 std::vector<std::vector<int32_t>> &permutations,
					 Covariates *cov,
					 TaskQueue &tq) {
  std::stringstream current;
  std::string gene;
  std::string transcript;
  std::string header;
  std::map<std::string, arma::uword> nvariants;
  int lineno = 0;
  int ngenes = 0;
  int ntranscripts = 0;

  std::ifstream ifs(tp.genotypes_path);
  std::string line;

  RJBUtil::Splitter<std::string> gene_list;
  // Parse gene list
  if (tp.gene_list) {
	gene_list = RJBUtil::Splitter<std::string>(*tp.gene_list, ",");
  }

  Bed bed;
  if (tp.bed)
	bed = Bed(*tp.bed);

  Weight casm;
  if (tp.weight)
	casm = Weight(*tp.weight);

  while (std::getline(ifs, line)) {
	if (lineno == 0) {
	  header = line;
	  if (!(*cov).is_sorted()) {
		cov->sort_covariates(header);
	  }

	  // Permute SKAT, SKATO, and BURDEN normally
	  Permute perm;
	  std::shared_ptr<std::vector<std::vector<int32_t>>> tmp_ptr(&permutations);
	  if (tp.stage_1_permutations > 0 && tp.method != "SKAT" && tp.method != "SKATO" && tp.method != "BURDEN") {
		perm.get_permutations(tmp_ptr,
							  cov->get_odds(),
							  cov->get_ncases(),
							  tp.stage_1_permutations,
							  tp.nthreads - 1);
	  }

	  // Permutation set output
	  if(tp.permute_set) {
		std::ofstream pset_ofs(*tp.permute_set);
		std::ofstream lr_ofs(*tp.permute_set + ".lr");
		for(const auto &p : permutations) {
		  for(const auto &v : p) {
		    pset_ofs << v << "\t";
		  }
		  pset_ofs << "\n";
		}
		pset_ofs.close();

		cov->get_probability().t().print(lr_ofs);
		lr_ofs.close();
		// Finish if there are no stage 2 permutations
		if(tp.stage_2_permutations == 0) {
		  std::exit(0);
		}
	  }

	  lineno++;
	  continue;
	}
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	if (splitter[0] == gene) {
	  // Add to current
	  RJBUtil::Splitter<std::string> var_split(splitter[2], "-");
	  if (tp.bed) {
		if (!bed.check_variant(var_split[0], var_split[1])) {
		  current << line << "\n";
		  if (splitter[1] == transcript) {
			nvariants[transcript]++;
		  } else {
			transcript = splitter[1];
			nvariants[transcript] = 1;
			ntranscripts++;
		  }
		}
	  } else {
		current << line << "\n";
		if (splitter[1] == transcript) {
		  nvariants[transcript]++;
		} else {
		  transcript = splitter[1];
		  nvariants[transcript] = 1;
		  ntranscripts++;
		}
	  }
	} else {
	  // Don't always create the gene
	  auto fit = std::find(gene_list.cbegin(), gene_list.cend(), gene);
	  bool skip_gene = (fit == gene_list.cend()) && tp.gene_list;
	  if (!gene.empty() && !skip_gene) {

		Gene gene_data(current, cov->get_nsamples(), nvariants, casm);
		Stage stage;

		// Choose starting stage
		if (tp.stage_1_permutations > 0) {
		  stage = Stage::Stage1;
		} else {
		  stage = Stage::Stage2;
		}

		// If we're looking at specific genes;
		if (tp.gene_list) {
		  if (gene_list.empty()) {
			// Gene list exhausted
			break;
		  }
		  // Found gene
		  if (fit != gene_list.cend()) {
			gene_list.erase(fit);
			// Ensure we have at least one variant for a submitted gene
			if (std::any_of(nvariants.cbegin(), nvariants.cend(), [&](const auto &v) { return v.second > 0; })) {
			  int total_s1_perm = 0;
			  int total_s2_perm = 0;
			  int total_success = 0;
			  for (int i = 0; i < tq.get_nthreads(); i++) {
				if (i == tq.get_nthreads() - 1) {
				  TaskArgs ta(stage,
							  gene_data,
							  *cov,
							  tp,
							  tp.success_threshold - total_success,
							  tp.stage_1_permutations - total_s1_perm,
							  tp.stage_2_permutations - total_s2_perm,
							  permutations);

				  // Limit adding jobs to prevent excessive memory usage
				  while (!tq.empty()) {
					std::this_thread::sleep_for(0.1s);
				  }
				  tq.dispatch(ta);
				} else {
				  TaskArgs ta(stage,
							  gene_data,
							  *cov,
							  tp,
							  tp.success_threshold / static_cast<int>(tq.get_nthreads()),
							  tp.stage_1_permutations / static_cast<int>(tq.get_nthreads()),
							  tp.stage_2_permutations / static_cast<int>(tq.get_nthreads()),
							  permutations);
				  // Add current permutations
				  total_s1_perm += tp.stage_1_permutations / tq.get_nthreads();
				  total_s2_perm += tp.stage_2_permutations / tq.get_nthreads();
				  total_success += tp.success_threshold / tq.get_nthreads();

				  // Limit adding jobs to prevent excessive memory usage
				  while (!tq.empty()) {
					std::this_thread::sleep_for(0.1s);
				  }
				  tq.dispatch(ta);
				}
			  }
			  ngenes++;
			} else {
			  std::cerr << "No variants in " << gene_data.get_gene() << ".\n";
			}
		  }
		} else {
		  // Ensure we have at least one variant for a submitted gene
		  if (std::any_of(nvariants.cbegin(), nvariants.cend(), [](const auto &v) { return v.second > 0; })) {
			TaskArgs ta(stage,
						gene_data,
						*cov,
						tp,
						permutations);
			// Limit adding jobs to prevent excessive memory usage
			while (!tq.empty()) {
			  std::this_thread::sleep_for(0.1s);
			}
			tq.dispatch(std::move(ta));
			ngenes++;
		  }
		}
	  }
	  if (tp.bed) {
		// Check if we add the line;
		RJBUtil::Splitter<std::string> var_split(splitter[2], "-");

		gene = splitter[0];
		transcript = splitter[1];

		// Reset current
		current.str("");
		current.clear();
		nvariants.clear();
		// Add header
		current << header << "\n";
		// New transcript
		ntranscripts++;

		if (!bed.check_variant(var_split[0], var_split[1])) {
		  current << line << "\n";
		  nvariants[transcript] = 1;
		} else {
		  // Exclude variant
		  nvariants[transcript] = 0;
		}
	  } else {
		gene = splitter[0];
		transcript = splitter[1];
		nvariants[transcript] = 1;
		ntranscripts++;

		current.str("");
		current.clear();
		current << header << "\n";
		current << line << "\n";
	  }
	}
  }
  auto fit = std::find(gene_list.cbegin(), gene_list.cend(), gene);
  bool skip_gene = (fit == gene_list.cend()) && tp.gene_list;

  // Submit final gene
  Gene gene_data(current, cov->get_nsamples(), nvariants, casm);
  Stage stage;

  // Choose starting stage
  if (tp.stage_1_permutations > 0) {
	stage = Stage::Stage1;
  } else {
	stage = Stage::Stage2;
  }
  // If we're looking at specific genes;
  if (tp.gene_list) {
	// Found gene
	if (std::find(gene_list.cbegin(), gene_list.cend(), gene) != gene_list.cend()) {
	  // Ensure we have at least one variant for a submitted gene
	  if (std::any_of(nvariants.cbegin(), nvariants.cend(), [&](const auto &v) { return v.second > 0; })) {
		// TODO Change logic for multiple jobs
		int total_s1_perm = 0;
		int total_s2_perm = 0;
		int total_success = 0;

		for (int i = 0; i < tq.get_nthreads(); i++) {
		  if (i == tq.get_nthreads() - 1) {
			TaskArgs ta(stage,
						gene_data,
						*cov,
						tp,
						tp.success_threshold - total_success,
						tp.stage_1_permutations - total_s1_perm,
						tp.stage_2_permutations - total_s2_perm,
						permutations);

			// Limit adding jobs to prevent excessive memory usage
			while (!tq.empty()) {
			  std::this_thread::sleep_for(0.1s);
			}
			tq.dispatch(ta);
		  } else {
			TaskArgs ta(stage,
						gene_data,
						*cov,
						tp,
						tp.success_threshold / static_cast<int>(tq.get_nthreads()),
						tp.stage_1_permutations / static_cast<int>(tq.get_nthreads()),
						tp.stage_2_permutations / static_cast<int>(tq.get_nthreads()),
						permutations);
			// Add current permutations
			total_s1_perm += tp.stage_1_permutations / tq.get_nthreads();
			total_s2_perm += tp.stage_2_permutations / tq.get_nthreads();
			total_success += tp.success_threshold / tq.get_nthreads();

			// Limit adding jobs to prevent excessive memory usage
			while (!tq.empty()) {
			  std::this_thread::sleep_for(0.1s);
			}
			tq.dispatch(ta);
		  }
		}
		ngenes++;
	  } else {
		std::cerr << "No variants in " << gene_data.get_gene() << ".\n";
	  }
	}
  } else {
	// Ensure at least one transcript
	if (std::any_of(nvariants.cbegin(), nvariants.cend(), [&](const auto &v) { return v.second > 0; })) {
	  TaskArgs ta(stage, gene_data, *cov, tp, permutations);

	  // Limit adding jobs to prevent excessive memory usage
	  while (!tq.empty()) {
		std::this_thread::sleep_for(0.1s);
	  }
	  tq.dispatch(std::move(ta));
	  ngenes++;
	}
  }

  // Wait for queue to finish processing
  tq.join();

  // Free permutation memory
  permutations.clear();

  if (tp.gene_list) {
	tq.duplicate();
	write_simple(tq, ntranscripts, ngenes, tp);
	write_detail(tq, tp);
  } else {
	write_simple(tq, ntranscripts, ngenes, tp);
	write_detail(tq, tp);
  }
}

