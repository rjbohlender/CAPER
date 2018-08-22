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
					 TaskQueue &tq) {
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

  Bed bed;
  if (!bed_file.empty())
	bed = Bed(bed_file);
  bool have_bed = !bed.empty();

  CASM casm;
  if (!casm_file.empty())
	casm = CASM(casm_file);

  while (std::getline(ifs, line)) {
	if (i==0) {
	  header = line;
	  i++;
	  continue;
	}
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	if (splitter[0]==gene) {
	  // Add to current
	  RJBUtil::Splitter<std::string> var_split(splitter[2], "-");
	  if (have_bed) {
		if (!bed.check_variant(var_split[0], var_split[1])) {
		  current << line << "\n";
		  if (splitter[1]==transcript) {
			nvariants[transcript]++;
		  } else {
			transcript = splitter[1];
			nvariants[transcript] = 1;
			ntranscripts++;
		  }
		}
	  } else {
		current << line << "\n";
		if (splitter[1]==transcript) {
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
					  permutations);
		  // Limit adding jobs to prevent excessive memory usage
		  while (!tq.empty()) {
			std::this_thread::sleep_for(0.1s);
		  }
		  tq.dispatch(std::move(ta));
		  ngenes++;
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

  // Ensure at least one transcript
  if (std::any_of(nvariants.cbegin(), nvariants.cend(), [&](const auto &v) { return v.second > 0; })) {
	TaskArgs
		ta(stage, gene_data, cov, successes, stage_1_perm, stage_2_perm, method, kernel, permutations);

	// Limit adding jobs to prevent excessive memory usage
	while (!tq.empty()) {
	  std::this_thread::sleep_for(0.1s);
	}
	tq.dispatch(std::move(ta));
	ngenes++;
  }

  // Wait for queue to finish processing
  tq.join();

  // Free permutation memory
  permutations.clear();

  // Print header and formatted results
  double permutation_mean = 0;
  double permutation_variance = 0;

  std::cout << std::setw(20) << "Gene";
  std::cout << std::setw(20) << "Transcript";
  std::cout << std::setw(20) << "Original";
  std::cout << std::setw(20) << "Empirical_P";
  std::cout << std::setw(20) << "Empirical_MidP";
  std::cout << std::setw(20) << "Successes";
  std::cout << std::setw(20) << "Permutations";
  std::cout << std::setw(20) << "MGIT";
  std::cout << std::setw(20) << "MGIT_Successes" << std::endl;
  for (auto &v : tq.get_results()) {
	std::cout << v.min_empirical_pvalue();

	for (const auto &k : v.get_gene().get_transcripts()) {
	  permutation_mean += v.results[k].permutations;
	}
  }
  permutation_mean /= ntranscripts;
  for (auto &v : tq.get_results()) {
	for (const auto &k : v.get_gene().get_transcripts()) {
	  permutation_variance += std::pow(v.results[k].permutations - permutation_mean, 2)/ntranscripts;
	}
  }
  std::cerr << "Permutation mean: " << permutation_mean << std::endl;
  std::cerr << "Permutation sd: " << std::sqrt(permutation_variance) << std::endl;
  std::cerr << "Genes submitted: " << ngenes << std::endl;
  std::cerr << "Transcripts submitted: " << ntranscripts << std::endl;
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
  parser.add_option("-n", "--nthreads")
	  .dest("nthreads")
	  .set_default(std::thread::hardware_concurrency()/2)
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
				  tq);
  return 0;
}