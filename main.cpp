#include "optparse.h"
#include "utility/split.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <thread>
#include <ctime>
#include <iomanip>
#include <algorithm>
#include <list>
#include <chrono>
#include <cassert>
#include <armadillo>

#include "data/covariates.hpp"
#include "data/taskqueue.hpp"
#include "third_party/stocc/stocc.h"
#include "third_party/stocc/randomc.h"
#include "statistics/methods.hpp"
#include "data/permutation.hpp"
#include "data/bed.hpp"
#include "data/casm.hpp"
#include "utility/filesystem.hpp"

using namespace std::chrono_literals;

void calc_mgit_pvalues(std::map<std::string, Result> &results, std::vector<std::string> &transcripts, const std::string &method) {
  // Skip MGIT if only a single transcript
  if (transcripts.size() == 1) {
    for(const auto &tr : transcripts) {
      auto &v = results[tr];
	  v.mgit_p = v.empirical_p;
	  v.mgit_successes = v.successes;
	}
	return;
  }

  // Shorthand
  unsigned long n = transcripts.size();
  int i, j, k;
  double successes;
  int max_perm = 0;

  // Get max_perm
  for(const auto &tr : transcripts) {
    if(results[tr].permutations > max_perm) {
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
	if (method == "SKATO") {
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

void print_results(TaskQueue &tq, int ntranscripts, int ngenes, bool genes) {
  std::cout << std::setw(20) << "Gene";
  std::cout << std::setw(20) << "Transcript";
  std::cout << std::setw(20) << "Original";
  std::cout << std::setw(20) << "Empirical_P";
  std::cout << std::setw(20) << "Empirical_MidP";
  std::cout << std::setw(20) << "Successes";
  std::cout << std::setw(20) << "Permutations";
  std::cout << std::setw(20) << "MGIT";
  std::cout << std::setw(20) << "MGIT_Successes" << std::endl;
  if(!genes) {
	// Print header and formatted results
	double permutation_mean = 0;
	double permutation_variance = 0;

	for (auto &v : tq.get_results()) {
	  std::cout << v.min_empirical_pvalue();

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

	  for(const auto &k : transcripts) {
		res[k] = std::move(v.results[k]);
	  }

	  // Check for additional results
	  if(tq.get_results().size() > 1) {
		std::vector<int> to_pop;
		for(int i = 1; i < tq.get_results().size(); i++) {
		  auto &cur_transcripts = tq.get_results()[i].get_gene().get_transcripts();
		  if(std::equal(cur_transcripts.cbegin(), cur_transcripts.cend(), transcripts.cbegin())) {
			// Mark for removal
			to_pop.push_back(i);
			// Combine
			for(const auto &k : cur_transcripts) {
			  res[k].combine(tq.get_results()[i].results[k]);
			}
		  }
		}
		// Remove remaining
		for(auto rit = to_pop.rbegin(); rit != to_pop.rend(); rit++) {
		  tq.get_results().erase(tq.get_results().begin() + *rit);
		}
		calc_mgit_pvalues(res, transcripts, v.get_methods().str());
	  }
	  tq.get_results().erase(tq.get_results().begin());
	  // Renew iterators
    } while(!tq.get_results().empty());

    for(const auto &v : res) {
      std::cout << v.second;
    }
  }
}

void initialize_jobs(const std::string &ifile,
					 const std::string &bed_file,
					 const std::string &casm_file,
					 const std::string &method,
					 const std::string &kernel,
					 size_t successes,
					 std::vector<std::vector<int32_t>> &permutations,
					 int stage_1_perm,
					 int stage_2_perm,
					 Covariates &cov,
					 TaskQueue &tq,
					 std::string gene_options,
					 bool genes,
					 bool adjust) {
  std::stringstream current;
  std::string gene;
  std::string transcript;
  std::string header;
  std::map<std::string, unsigned long> nvariants;
  int i = 0;
  int ngenes = 0;
  int ntranscripts = 0;

  std::ifstream ifs(ifile);
  std::string line;

  RJBUtil::Splitter<std::string> gene_list;
  // Parse gene list
  if (genes) {
	gene_list = RJBUtil::Splitter<std::string>(gene_options, ",");
  }

  Bed bed;
  if (!bed_file.empty())
	bed = Bed(bed_file);
  bool have_bed = !bed.empty();

  CASM casm;
  if (!casm_file.empty())
	casm = CASM(casm_file);

  while (std::getline(ifs, line)) {
	if (i == 0) {
	  header = line;
	  i++;
	  continue;
	}
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	if (splitter[0] == gene) {
	  // Add to current
	  RJBUtil::Splitter<std::string> var_split(splitter[2], "-");
	  if (have_bed) {
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
	  if (!gene.empty()) {
		Gene gene_data(current, cov.get_nsamples(), nvariants, casm);
		Stage stage;

		// Choose starting stage
		if (stage_1_perm > 0) {
		  stage = Stage::Stage1;
		} else {
		  stage = Stage::Stage2;
		}

		// If we're looking at specific genes;
		if (genes) {
		  // Found gene
		  if (std::find(gene_list.cbegin(), gene_list.cend(), gene_data.get_gene()) != gene_list.cend()) {
			// Ensure we have at least one variant for a submitted gene
			if (std::any_of(nvariants.cbegin(), nvariants.cend(), [&](const auto &v) { return v.second > 0; })) {
			  // TODO Change logic for multiple jobs
			  int total_s1_perm = 0;
			  int total_s2_perm = 0;
			  for (int i = 0; i < tq.get_nthreads(); i++) {
				if (i == tq.get_nthreads() - 1) {
				  TaskArgs ta(stage,
							  gene_data,
							  cov,
							  successes,
							  stage_1_perm - total_s1_perm,
							  stage_2_perm - total_s2_perm,
							  method,
							  kernel,
							  permutations,
							  adjust);

				  // Limit adding jobs to prevent excessive memory usage
				  while (!tq.empty()) {
					std::this_thread::sleep_for(0.1s);
				  }
				  tq.dispatch(ta);
				} else {
				  TaskArgs ta(stage,
							  gene_data,
							  cov,
							  successes,
							  stage_1_perm / tq.get_nthreads(),
							  stage_2_perm / tq.get_nthreads(),
							  method,
							  kernel,
							  permutations,
							  adjust);
				  // Add current permutations
				  total_s1_perm += stage_1_perm / tq.get_nthreads();
				  total_s2_perm += stage_2_perm / tq.get_nthreads();

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
		  if (std::any_of(nvariants.cbegin(), nvariants.cend(), [&](const auto &v) { return v.second > 0; })) {
			TaskArgs ta(stage,
						gene_data,
						cov,
						successes,
						stage_1_perm,
						stage_2_perm,
						method,
						kernel,
						permutations,
						adjust);
			// Limit adding jobs to prevent excessive memory usage
			while (!tq.empty()) {
			  std::this_thread::sleep_for(0.1s);
			}
			tq.dispatch(std::move(ta));
			ngenes++;
		  }
		}
	  }
	  if (have_bed) {
		// Check if we add the line;
		RJBUtil::Splitter<std::string> var_split(splitter[2], "-");

		gene = splitter[0];
		transcript = splitter[1];

		// Reset current
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

		current.clear();
		current << header << "\n";
		current << line << "\n";
	  }
	}
  }

  // Submit final gene
  Gene gene_data(current, cov.get_nsamples(), nvariants, casm);
  Stage stage;

  // Choose starting stage
  if (stage_1_perm > 0) {
	stage = Stage::Stage1;
  } else {
	stage = Stage::Stage2;
  }
  // If we're looking at specific genes;
  if (genes) {
	// Found gene
	if (std::find(gene_list.cbegin(), gene_list.cend(), gene_data.get_gene()) != gene_list.cend()) {
	  // Ensure we have at least one variant for a submitted gene
	  if (std::any_of(nvariants.cbegin(), nvariants.cend(), [&](const auto &v) { return v.second > 0; })) {
		// TODO Change logic for multiple jobs
		int total_s1_perm = 0;
		int total_s2_perm = 0;
		for (int i = 0; i < tq.get_nthreads(); i++) {
		  if (i == tq.get_nthreads() - 1) {
			TaskArgs ta(stage,
						gene_data,
						cov,
						successes,
						stage_1_perm - total_s1_perm,
						stage_2_perm - total_s2_perm,
						method,
						kernel,
						permutations,
						adjust);

			// Limit adding jobs to prevent excessive memory usage
			while (!tq.empty()) {
			  std::this_thread::sleep_for(0.1s);
			}
			tq.dispatch(ta);
		  } else {
			TaskArgs ta(stage,
						gene_data,
						cov,
						successes,
						stage_1_perm / tq.get_nthreads(),
						stage_2_perm / tq.get_nthreads(),
						method,
						kernel,
						permutations,
						adjust);
			// Add current permutations
			total_s1_perm += stage_1_perm / tq.get_nthreads();
			total_s2_perm += stage_2_perm / tq.get_nthreads();

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
	  TaskArgs
		  ta(stage,
			 gene_data,
			 cov,
			 successes,
			 stage_1_perm,
			 stage_2_perm,
			 method,
			 kernel,
			 permutations,
			 adjust);

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

  print_results(tq, ntranscripts, ngenes, genes);
}


int main(int argc, char **argv) {
  // Only using C++ IO.
  std::ios_base::sync_with_stdio(false);

  std::list<std::string> method_choices = {
	  "CALPHA",
	  "CMC",
	  "VT",
	  "WSS",
	  "SKAT",
	  "SKATO",
	  "VAAST"
  };

  optparse::OptionParser parser =
	  optparse::OptionParser().description("Permutation tool for gene-based analysis.");

  parser.add_option("-g", "--genotypes")
	  .dest("genotypes")
	  .nargs(1)
	  .help("The genotype table file.");
  parser.add_option("-c", "--covariates")
	  .dest("covariates")
	  .nargs(1)
	  .help("The covariate table file, including phenotypes.\n"
			"Phenotypes {0=control, 1=case} should be in the first column.");
  parser.add_option("-b", "--bed_file")
	  .dest("bed_file")
	  .set_default("")
	  .help("A bed file to be used as a filter. All specified regions will be excluded.");
  parser.add_option("-w", "--casm_file")
	  .dest("casm_file")
	  .set_default("")
	  .help("A file providing weights for VAAST.");
  parser.add_option("-t", "--nthreads")
	  .dest("nthreads")
	  .set_default(std::thread::hardware_concurrency() / 2)
	  .type("size_t")
	  .help("The number of threads. The default is half of the cpu count.");
  parser.add_option("-1", "--stage_1_max_perm")
	  .dest("stage_1_max_perm")
	  .set_default(10000)
	  .help("The maximum number of permutations to be performed in the first stage. The default is 10,000.");
  parser.add_option("-2", "--stage_2_max_perm")
	  .dest("stage_2_max_perm")
	  .set_default(100000)
	  .help("The maximum number of permutations to be performed in the second stage.\n"
			"The default is 100,000. This includes the number performed in the first stage.");
  parser.add_option("-m", "--method")
	  .dest("method")
	  .set_default("VAAST")
	  .choices(method_choices.begin(), method_choices.end())
	  .help("The statistical method to be used.\nOptions: {CALPHA, CMC, SKAT, WSS, VAAST, VT}.\nThe default is VAAST.");
  parser.add_option("-n", "--no_adjust")
	  .dest("adjust")
	  .action("store_false")
	  .set_default("1")
	  .help("Disable small sample size adjustment for SKATO.");
  parser.add_option("-q", "--quiet")
	  .dest("verbose")
	  .action("store_false")
	  .set_default("1")
	  .help("Don't print status messages.");
  parser.add_option("-s", "--successes")
	  .dest("successes")
	  .set_default(200)
	  .type("size_t")
	  .help("Number of successes for early termination.");
  parser.add_option("-k", "--kernel")
	  .dest("kernel")
	  .set_default("Linear")
	  .help("Kernel for use with SKAT.\nOne of: {Linear, wLinear, IBS, wIBS, Quadratic, twoWayX}.\nDefault is Linear.");
  parser.add_option("-l", "--genes")
	  .dest("genes")
	  .set_default("")
	  .help("A comma separated list of genes to analyze.");
  parser.add_help_option();

  const optparse::Values options = parser.parse_args(argc, argv);

  // Print help and quit if missing input files
  if (!options.is_set("genotypes") || !options.is_set("covariates")) {
	parser.print_help();
	std::exit(1);
  }

  size_t nthreads = options.get("nthreads");
  assert(nthreads > 1);

  if (options.get("verbose")) {
	std::cerr << "genotypes: " << options["genotypes"] << "\n";
	std::cerr << "covariates: " << options["covariates"] << "\n";
	std::cerr << "bed_file: " << options["bed_file"] << "\n";
	std::cerr << "casm_file: " << options["casm_file"] << "\n";
	std::cerr << "success threshold: " << options["successes"] << "\n";
	std::cerr << "nthreads: " << options["nthreads"] << "\n";
	std::cerr << "stage_1_max_perm: " << options["stage_1_max_perm"] << "\n";
	std::cerr << "stage_2_max_perm: " << options["stage_2_max_perm"] << "\n";
  }
  TaskQueue tq(nthreads - 1, options.get("verbose"));

  // Check for correct file paths
  if (!check_file_exists(options["genotypes"])) {
	std::cerr << "Incorrect file path for genotypes." << std::endl;
	parser.print_help();
	std::exit(1);
  }
  if (!check_file_exists(options["covariates"])) {
	std::cerr << "Incorrect file path for covariates." << std::endl;
	parser.print_help();
	std::exit(1);
  }
  if (!options["bed_file"].empty() && !check_file_exists(options["bed_file"])) {
	std::cerr << "Incorrect file path for bed_file." << std::endl;
	parser.print_help();
	std::exit(1);
  }
  if (!options["casm_file"].empty() && !check_file_exists(options["casm_file"])) {
	std::cerr << "Incorrect file path for casm_file." << std::endl;
	parser.print_help();
	std::exit(1);
  }

  Covariates cov(options["covariates"]);
  std::vector<std::vector<int32_t>> permutations;
  int s1_perm = options.get("stage_1_max_perm");
  Permute perm;
  // Permute SKAT and SKATO normally
  if (s1_perm > 0 && options["method"] != "SKAT" && options["method"] != "SKATO") {
	permutations = perm.get_permutations(options.get("stage_1_max_perm"), cov.get_odds(), cov.get_ncases());
  }

  // Initialize randomization for SKAT
  //if(options["method"] == "SKAT" || options["method"] == "SKATO")
  arma::arma_rng::set_seed_random();

  initialize_jobs(options["genotypes"],
				  options["bed_file"],
				  options["casm_file"],
				  options["method"],
				  options["kernel"],
				  options.get("successes"),
				  permutations,
				  options.get("stage_1_max_perm"),
				  options.get("stage_2_max_perm"),
				  cov,
				  tq,
				  options["genes"],
				  options.is_set_by_user("genes"),
				  options.get("adjust"));
  return 0;
}