#include "../utility/split.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <thread>
#include <ctime>
#include <iomanip>
#include <set>
#include <cassert>

#include <boost/program_options.hpp>
#include <boost/optional.hpp>

#include "stocc/stocc.h"
#include "stocc/randomc.h"
#include "../statistics/methods.hpp"
#include "../data/permutation.hpp"
#include "../utility/filesystem.hpp"
#include "../utility/main_support.hpp"
#include "../utility/jobdispatcher.hpp"
#include "../utility/taskparams.hpp"
#include "../utility/reporter.hpp"
#include "../power/powerop.hpp"

namespace po = boost::program_options;

int main(int argc, char **argv) {
  // Only using C++ IO.
  std::ios_base::sync_with_stdio(false);

  // Run timer
  arma::wall_clock timer;
  timer.tic();

  // BOOST Program Options Implementation
  po::options_description all("");
  po::options_description visible("Permutation tool for gene-based rare-variant analysis.\nAllowed options");
  po::options_description required("Required");
  po::options_description optional("Optional");
  po::options_description vaast("VAAST Options");
  po::options_description skat("SKAT Options");
  po::options_description cmc("CMC Options");
  po::options_description hidden("Hidden Options");
  po::variables_map vm;

  bool verbose = true;
  bool adjust = true;
  bool score_only_minor = false;
  bool score_only_alternative = false;
  bool testable = false;
  bool linear = false;
  bool nodetail = false;
  bool top_only = false;
  bool biallelic = false;
  bool covadj = true;
  std::vector<int> gene_range;
  std::vector<std::string> power;
  boost::optional<std::string> bed;
  boost::optional<std::string> weight;
  boost::optional<std::string> gene_list;
  boost::optional<std::string> permute_set;
  boost::optional<double> pthresh;
  boost::optional<arma::uword> approximate;

  try {
    required.add_options()
				("input,i",
				 po::value<std::string>()->required(),
				 "Genotype matrix file path.")
				("covariates,c",
				 po::value<std::string>()->required(),
				 "The covariate matrix file, tab separated.\nFormat = sample_id cov1 ...")
				("ped,p",
				 po::value<std::string>()->required(),
				 "Path to the .ped file containing the sample phenotypes.")
				("output,o",
				 po::value<std::string>()->required(),
				 "Path to output directory. Two files will be output: a simple transcript level results file, and a detailed variant level result file.")
				 ;
    optional.add_options()
				("bed-file,b",
				 po::value(&bed),
				 "A bed file to be used as a filter. All specified regions will be excluded.")
				("weight-file,w",
				 po::value(&weight),
				 "A file providing weights.")
				("nthreads,t",
				 po::value<size_t>()->default_value(std::thread::hardware_concurrency() / 2 + 1),
				 "The number of threads. Minimum number of threads = 2. n + 1 threads, with one parent thread and n worker threads.")
				("method,m",
				 po::value<std::string>()->default_value("VAAST"),
				 "The statistical method to be used.\nOptions: {BURDEN, CALPHA, CMC, SKAT, WSS, VAAST, VT}.")
				("range",
				 po::value(&gene_range)->multitoken(),
				 "A range of genes to analyze from the matrix file. Takes two values, a start gene number, and end gene number.\nThe program will only provide results for the values in that range. Gene count starts at 1.\nUseful for starting multiple jobs on a cluster each processing part of a file.")
				("stage_1_max_perm,1",
				 po::value<arma::uword>()->default_value(0),
				 "The maximum number of permutations to be performed in the first stage permutation. A small number is recommended if your sample is large.")
				("stage_2_max_perm,2",
				 po::value<arma::uword>()->default_value(1000000),
				 "The maximum number of permutations to be performed in the second stage, the collapsing step.")
				("mac",
				 po::value<arma::uword>()->default_value(250),
				 "Minor allele count cutoff.")
				("maf,r",
				 po::value<double>()->default_value(0.05),
				 "Minor allele frequency cutoff. Default equivalent to XQC.")
				("pthresh,j",
				 po::value(&pthresh),
				 "The threshold to terminate permutation based on whether it is outside the p-value CI.")
				("top_only",
				 po::bool_switch(&top_only),
				 "Output only the top transcript in the simple file.")
				("successes,s",
				 po::value<arma::uword>()->default_value(200),
				 "Number of successes for early termination.")
				("genes,l",
				 po::value(&gene_list),
				 "A comma-separated list of genes to analyze.")
				("nodetail",
				 po::bool_switch(&nodetail),
				 "Don't produce detailed, variant level output.")
				("approx,a",
				 po::value(&approximate),
				 "Group minor allele carriers into n bins with shared average odds. This option is useful for very large data sets where the total number of minor allele carriers to be permuted can be very large, and result in extreme run times.")
				 ;
    vaast.add_options()
			 ("group_size,g",
			  po::value<arma::uword>()->default_value(0),
			  "Group size. VAAST can collapse variants into groups of variants with adjacent weights.")
			 ("score_only_minor",
			  po::bool_switch(&score_only_minor),
			  "Score only minor alleles in VAAST.")
			 ("score_only_alternative",
			  po::bool_switch(&score_only_alternative),
			  "Score only alternative alleles in VAAST.")
			 ("testable",
			  po::bool_switch(&testable),
			  "Return scores only for genes with at least scoreable variants in VAAST.")
			 ("biallelic",
			  po::bool_switch(&biallelic),
			  "Additional term for biallelic variants. For detecting potentially recessive variants.")
			  ;
    skat.add_options()
			("kernel,k",
			 po::value<std::string>()->default_value("wLinear"),
			 "Kernel for use with SKAT.\nOne of: {Linear, wLinear}.")
			("qtl",
			 po::bool_switch(&linear),
			 "Analyze a quantitative trait. Values are assumed to be finite floating point values.")
			("beta_weights",
			 po::value<std::string>()->default_value("1,25"),
			 "Parameters for the beta distribution. Two values, comma separated corresponding to a,b.")
			 ;
    cmc.add_options()
		   ("cmcmaf",
		   	po::value<double>()->default_value(0.005),
		   	"Minor allele frequency cutoff for CMC collapsing.")
		   ;
	all.add_options()
			("help,h", "Print this help message.")
			("quiet,q", "Don't print status messages.")
			;
	hidden.add_options()
			// ("no_adjust,n", "Disable small sample size adjustment for SKATO.")
			("permute_out",
			 po::value(&permute_set),
			 "Output permutations to the given file.")
			("nocovadj",
			 po::bool_switch(&covadj),
			 "Do Fisher-Yates shuffling instead of covariate adjusted permutation.")
			;
	all.add(required).add(optional).add(vaast).add(skat).add(cmc).add(hidden);
	visible.add(required).add(optional).add(vaast).add(skat).add(cmc);
	po::store(po::parse_command_line(argc, argv, all), vm);
	if (vm.count("help")) {
	  std::cerr << visible << "\n";
	  return 1;
	}

	std::vector<int> range_opt;
	if (!vm["range"].empty() && ((range_opt = vm["range"].as<std::vector<int>>()).size() != 2 || (!vm["range"].empty() && !vm["genes"].empty()))) {
	  std::cerr << "--range takes two integer arguments\n";
	  std::cerr << "--range cannot be used with the --genes or -l option.\n";
	  std::cerr << visible << "\n";
	  return 1;
	}

	po::notify(vm);

	if (vm.count("quiet")) {
	  verbose = false;
	}
	if (vm.count("no_adjust")) {
	  adjust = false;
	}
  } catch (po::required_option &e) {
	std::cerr << "Missing required option:\n" << e.what() << "\n";
	std::cerr << visible << "\n";
	return 1;
  } catch (std::exception &e) {
	std::cerr << "Error: " << e.what() << "\n";
	return 1;
  }

  std::set<std::string> method_choices = {
	  "BURDEN",
	  "CALPHA",
	  "CMC",
	  "VT",
	  "WSS",
	  "RVT1",
	  "RVT2",
	  "SKAT",
	  "SKATO",
	  "VAAST"
  };

  std::set<std::string> kernel_choices = {
	  "Linear",
	  "wLinear"
	  //"IBS",
	  //"wIBS",
	  //"Quadratic",
	  //"twoWayX"
  };

  if (method_choices.count(vm["method"].as<std::string>()) == 0) {
	// Method not among choices
	std::cerr << "Method must be one of {BURDEN, CALPHA, CMC, RVT1, RVT2, SKAT, SKATO, VAAST, VT, WSS}.\n";
	std::cerr << visible << "\n";
	return 1;
  }

  if (kernel_choices.count(vm["kernel"].as<std::string>()) == 0) {
	// Method not among choices
	std::cerr << "Kernel must be one of {Linear, wLinear}.\n";
	std::cerr << visible << "\n";
	return 1;
  }

  /**********************
   * Setup task parameters
   **********************/
  TaskParams tp;

  RJBUtil::Splitter<std::string> beta_split(vm["beta_weights"].as<std::string>(), ",");

  // Store full command
  std::stringstream cmd_ss;
  for (int i = 0; i < argc; i++) {
	if (i == argc - 1) {
	  cmd_ss << argv[i];
	} else {
	  cmd_ss << argv[i] << " ";
	}
  }

  tp.base = argv[0];

  tp.full_command = cmd_ss.str();

  tp.success_threshold = vm["successes"].as<arma::uword>();
  tp.stage_1_permutations = vm["stage_1_max_perm"].as<arma::uword>();
  tp.stage_2_permutations = vm["stage_2_max_perm"].as<arma::uword>();
  tp.total_permutations = std::max(tp.stage_1_permutations, tp.stage_2_permutations);
  tp.method = vm["method"].as<std::string>();
  // File paths and option status
  tp.program_path = argv[0];
  tp.genotypes_path = vm["input"].as<std::string>();
  tp.covariates_path = vm["covariates"].as<std::string>();
  tp.ped_path = vm["ped"].as<std::string>();
  tp.output_path = vm["output"].as<std::string>();
  tp.maf = vm["maf"].as<double>();
  tp.cmcmaf = vm["cmcmaf"].as<double>();
  tp.group_size = vm["group_size"].as<arma::uword>();
  tp.score_only_minor = score_only_minor;
  tp.score_only_alternative = score_only_alternative;
  tp.bed = bed;
  tp.weight = weight;
  tp.permute_set = permute_set;
  // Threading
  tp.nthreads = vm["nthreads"].as<size_t>();
  // Options
  tp.verbose = verbose;
  tp.gene_list = gene_list;
  tp.nodetail = nodetail;
  tp.top_only = top_only;
  tp.mac = vm["mac"].as<arma::uword>();
  tp.pthresh = pthresh;
  tp.approximate = approximate;
  // SKAT Options
  tp.kernel = vm["kernel"].as<std::string>();
  tp.adjust = adjust;
  tp.linear = linear;
  // Beta weights
  tp.a = std::stoi(beta_split[0]);
  tp.b = std::stoi(beta_split[1]);
  // Testability
  tp.testable = testable;
  tp.biallelic = biallelic;

  tp.covadj = covadj;

  tp.alternate_permutation = tp.method == "SKATO" || tp.method == "SKAT" || tp.method == "BURDEN" || tp.method == "VT";
  tp.quantitative = tp.method == "RVT1" || tp.method == "RVT2" || tp.method == "SKATO" || tp.method == "SKAT" || tp.method == "BURDEN" || tp.method == "VT";
  tp.analytic = tp.method == "SKATO" || (tp.method == "SKAT" && tp.total_permutations == 0) || tp.method == "RVT1" || tp.method == "RVT2" || tp.method == "CMC";
  if(tp.linear && !tp.quantitative) {
    std::cerr << "Quantitative trait analysis is only supported for the RVT1, RVT2, SKATO, SKAT, and BURDEN methods." << std::endl;
    std::exit(1);
  }

  std::vector<int> range_opt;
  if(!vm["range"].empty() && (range_opt = vm["range"].as<std::vector<int>>()).size() == 2) {
    tp.range_start = range_opt[0];
    tp.range_end = range_opt[1];
  }

  if(tp.mac <= 0) {
    std::cerr << "Minor allele count cutoff must be greater than zero." << std::endl;
    std::exit(1);
  } else if(tp.mac > 500) {
    std::cerr << "WARNING: This software is concerned with evaluating rare events. With a minor allele cutoff > 500, you should consider analyzing those variants using single marker tests." << std::endl;
  }

  if (tp.nthreads < 2) {
    std::cerr << "Thread count must be >= 2." << std::endl;
    std::cerr << visible << std::endl;
  }
  if(tp.permute_set && tp.nthreads) {
    std::cerr << "Restricting to a single worker thread. Permute set output is not threadsafe." << std::endl;
    tp.nthreads = 2;
  }

  if (tp.method == "SKATO") {
	// Different defaults for SKATO
	if(vm["stage_1_max_perm"].defaulted() || vm["stage_1_max_perm"].as<arma::uword>() > 0) {
	  tp.stage_1_permutations = 0;
	}
	if(vm["stage_2_max_perm"].defaulted()) {
	  tp.stage_2_permutations = 0;
	}
  }

  if (tp.verbose) {
	std::cerr << "genotypes: " << tp.genotypes_path << "\n";
	std::cerr << "covariates: " << tp.covariates_path << "\n";
	std::cerr << "ped: " << tp.ped_path << "\n";
	if(tp.bed)
	  std::cerr << "bed_file: " << *tp.bed << "\n";
	if(tp.weight)
	  std::cerr << "weight_file: " << *tp.weight << "\n";
	std::cerr << "output: " << tp.output_path << "\n";
	std::cerr << "method: " << tp.method << "\n";
	std::cerr << "success threshold: " << tp.success_threshold << "\n";
	std::cerr << "nthreads: " << tp.nthreads << "\n";
	std::cerr << "stage_1_max_perm: " << tp.stage_1_permutations << "\n";
	std::cerr << "stage_2_max_perm: " << tp.stage_2_permutations << "\n";
  }

  // Check for correct file paths
  if (!check_file_exists(tp.genotypes_path)) {
	std::cerr << "Incorrect file path for genotypes." << std::endl;
	std::cerr << visible << "\n";
	std::exit(1);
  }
  if (!check_file_exists(tp.covariates_path)) {
	std::cerr << "Incorrect file path for covariates." << std::endl;
	std::cerr << visible << "\n";
	std::exit(1);
  }
  if (!check_file_exists(tp.ped_path)) {
	std::cerr << "Incorrect file path for ped." << std::endl;
	std::cerr << visible << "\n";
	std::exit(1);
  }
  if (tp.bed && !check_file_exists(*tp.bed)) {
	std::cerr << "Incorrect file path for bed_file." << std::endl;
	std::cerr << visible << "\n";
	std::exit(1);
  }
  if (tp.weight && !check_file_exists(*tp.weight)) {
	std::cerr << "Incorrect file path for weight_file." << std::endl;
	std::cerr << visible << "\n";
	std::exit(1);
  }
  if (!check_directory_exists(tp.output_path)) {
  	std::cerr << "Output path is invalid." << std::endl;
  	std::cerr << visible << "\n";
  	std::exit(1);
  }
  // Initialize randomization
  arma::arma_rng::set_seed_random();

  std::shared_ptr<Reporter> reporter = nullptr;
  reporter = std::make_shared<Reporter>(tp);
  JobDispatcher<CARVAOp, CARVATask, Reporter> jd(tp, reporter);

  double n = timer.toc();
  std::cerr << "Elapsed time: " << n << std::endl;
  return 0;
}