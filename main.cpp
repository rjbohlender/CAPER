#include "utility/split.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <thread>
#include <ctime>
#include <iomanip>
#include <algorithm>
#include <set>
#include <cassert>

#include <boost/program_options.hpp>
#include <boost/optional.hpp>

#include "third_party/stocc/stocc.h"
#include "third_party/stocc/randomc.h"
#include "statistics/methods.hpp"
#include "data/permutation.hpp"
#include "utility/filesystem.hpp"
#include "utility/main_support.hpp"

namespace po = boost::program_options;

int main(int argc, char **argv) {
  // Only using C++ IO.
  std::ios_base::sync_with_stdio(false);

  // BOOST Program Options Implementation
  po::options_description desc("Permutation tool for gene-based rare-variant analysis.\nAllowed options");
  po::variables_map vm;

  bool verbose = true;
  bool adjust = true;
  boost::optional<std::string> bed;
  boost::optional<std::string> casm;
  boost::optional<std::string> genes;

  try {
	desc.add_options()
			("help,h", "Print this help message.")
			("genotypes,g", po::value<std::string>()->required(), "Genotype matrix file.")
			("covariates,c", po::value<std::string>()->required(), "The covariate table file, including phenotypes.\nPhenotypes {0=control, 1=case} should be in the first column.")
			("bed-file,b", po::value(&bed), "A bed file to be used as a filter. All specified regions will be excluded.")
			("casm-file,w", po::value(&casm), "A file providing weights for VAAST.")
			("nthreads,t",
			 po::value<size_t>()->default_value(std::thread::hardware_concurrency() / 2),
			 "The number of threads. The default is half of the cpu count.")
			("method,m", po::value<std::string>()->default_value("VAAST"), "The statistical method to be used.\nOptions: {CALPHA, CMC, SKAT, WSS, VAAST, VT}.\nThe default is VAAST.")
			("stage_1_max_perm,1", po::value<int>()->default_value(100000), "The maximum number of permutations to be performed in the first stage. The default is 100,000.")
			("stage_2_max_perm,2", po::value<int>()->default_value(1000000), "The maximum number of permutations to be performed in the second stage. The default is 1,000,000.")
			("no_adjust,n", "Disable small sample size adjustment for SKATO.")
			("kernel,k", po::value<std::string>()->default_value("wLinear"), "Kernel for use with SKAT.\nOne of: {Linear, wLinear, IBS, wIBS, Quadratic, twoWayX}. Default is wLinear.")
			("beta_weights", po::value<std::string>()->default_value("1,25"), "Parameters for the beta distribution. Two values, comma separated corresponding to a,b. Default is 1,25.")
			("successes,s", po::value<int>()->default_value(200), "Number of successes for early termination.")
			("genes,l", po::value(&genes), "A comma-separated list of genes to analyze.")
			("quiet,q", "Don't print status messages.")
		;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if(vm.count("help")) {
	  std::cerr << desc << "\n";
	  return 1;
	}

	po::notify(vm);

	if(vm.count("quiet")) {
	  verbose = false;
	}
	if(vm.count("no_adjust")) {
	  adjust = false;
	}
  } catch(po::required_option &e) {
    std::cerr << "Missing required option:\n" << e.what() << "\n";
    std::cerr << desc << "\n";
    return 1;
  } catch(std::exception &e) {
	std::cerr << "Error: " << e.what() << "\n";
	return 1;
  }

  std::set<std::string> method_choices = {
	  "CALPHA",
	  "CMC",
	  "VT",
	  "WSS",
	  "SKAT",
	  "SKATO",
	  "VAAST"
  };

  std::set<std::string> kernel_choices = {
  	"Linear",
  	"wLinear",
  	"IBS",
  	"wIBS",
  	"Quadratic",
  	"twoWayX"
  };

  if(method_choices.count(vm["method"].as<std::string>()) == 0) {
    // Method not among choices
    std::cerr << "Method must be one of {CALPHA, CMC, VT, WSS, SKAT, SKATO, VAAST}.\n";
    std::cerr << desc << "\n";
    return 1;
  }

  if(kernel_choices.count(vm["kernel"].as<std::string>()) == 0) {
	// Method not among choices
	std::cerr << "Kernel must be one of {Linear, wLinear, IBS, wIBS, Quadratic, twoWayX}.\n";
	std::cerr << desc << "\n";
	return 1;
  }

  // Setup task parameters
  TaskParams tp{};

  RJBUtil::Splitter<std::string> beta_split(vm["beta_weights"].as<std::string>(), ",");

  tp.success_threshold = vm["successes"].as<int>();
  tp.stage_1_permutations = vm["stage_1_max_perm"].as<int>();
  tp.stage_2_permutations = vm["stage_2_max_perm"].as<int>();
  tp.total_permutations = std::max(tp.stage_1_permutations, tp.stage_2_permutations);
  tp.method = vm["method"].as<std::string>();
  // File paths and option status
  tp.program_path = argv[0];
  tp.genotypes_path = vm["genotypes"].as<std::string>();
  tp.covariates_path = vm["covariates"].as<std::string>();
  if(bed) {
	tp.bed = true;
	tp.bed_path = *bed;
  } else {
    tp.bed = false;
    tp.bed_path = "";
  }
  if(casm) {
    tp.casm = true;
    tp.casm_path = *casm;
  } else {
    tp.casm = false;
    tp.casm_path = "";
  }
  // Threading
  tp.nthreads = vm["nthreads"].as<size_t>();
  // Options
  tp.verbose = verbose;
  if(genes) {
	tp.genes = true;
	tp.gene_list = *genes;
  } else {
    tp.genes = false;
    tp.gene_list = "";
  }
  // SKAT Options
  tp.kernel = vm["kernel"].as<std::string>();
  tp.adjust = adjust;
  // Beta weights
  tp.a = std::stod(beta_split[0]);
  tp.b = std::stod(beta_split[1]);

  assert(tp.nthreads > 1);

  if (tp.verbose) {
	std::cerr << "genotypes: " << tp.genotypes_path << "\n";
	std::cerr << "covariates: " << tp.covariates_path << "\n";
	std::cerr << "bed_file: " << tp.bed_path << "\n";
	std::cerr << "casm_file: " << tp.casm_path << "\n";
	std::cerr << "success threshold: " << tp.success_threshold << "\n";
	std::cerr << "nthreads: " << tp.nthreads << "\n";
	std::cerr << "stage_1_max_perm: " << tp.stage_1_permutations << "\n";
	std::cerr << "stage_2_max_perm: " << tp.stage_2_permutations << "\n";
	if (tp.method == "SKATO") {
	  std::cerr << "adjust: " << tp.adjust << "\n";
	}
  }

  TaskQueue tq(tp.nthreads - 1, tp.verbose);

  // Check for correct file paths
  if (!check_file_exists(tp.genotypes_path)) {
	std::cerr << "Incorrect file path for genotypes." << std::endl;
	std::cerr << desc << "\n";
	std::exit(1);
  }
  if (!check_file_exists(tp.covariates_path)) {
	std::cerr << "Incorrect file path for covariates." << std::endl;
	std::cerr << desc << "\n";
	std::exit(1);
  }
  if (tp.bed && !check_file_exists(tp.bed_path)) {
	std::cerr << "Incorrect file path for bed_file." << std::endl;
	std::cerr << desc << "\n";
	std::exit(1);
  }
  if (tp.casm && !check_file_exists(tp.casm_path)) {
	std::cerr << "Incorrect file path for casm_file." << std::endl;
	std::cerr << desc << "\n";
	std::exit(1);
  }

  Covariates cov(tp.covariates_path);
  std::vector<std::vector<int32_t>> permutations;

  Permute perm;
  // Permute SKAT and SKATO normally
  if (tp.stage_1_permutations > 0 && tp.method != "SKAT" && tp.method != "SKATO") {
	permutations = perm.get_permutations(tp.stage_1_permutations, cov.get_odds(), cov.get_ncases());
  }

  // Initialize randomization
  arma::arma_rng::set_seed_random();

  initialize_jobs(tp, permutations, cov, tq);
  return 0;
}