#include "../utility/split.hpp"

#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <thread>
#include <set>
#include <cassert>

#include <boost/program_options.hpp>
#include <boost/optional.hpp>
#include <unistd.h>

#include "../statistics/methods.hpp"
#include "../utility/filesystem.hpp"
#include "../utility/jobdispatcher.hpp"
#include "../power/powerop.hpp"

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

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
  bool linear = false;
  bool no_detail = false;
  bool top_only = false;
  bool biallelic = false;
  bool nocovadj = false;
  bool stats = false;
  bool alternate_grouping = false;
  bool no_weights = false;
  bool impute_to_mean = false;
  bool ma_count = true;
  bool whole_gene = false;
  bool saddlepoint = false;
  bool hotellings = false;
  bool wald = false;
  bool var_collapsing = false;
  std::vector<int> gene_range;
  std::vector<std::string> power;
  boost::optional<std::string> bed;
  boost::optional<std::string> weights;
  boost::optional<std::string> gene_list;
  boost::optional<std::string> feature_list;
  boost::optional<std::string> permute_set;
  boost::optional<double> pthresh;
  boost::optional<arma::uword> approximate;
  boost::optional<int> seed;
  boost::optional<arma::uword> max_perms;
  boost::optional<double> testable;
  double soft_maf_filter = 1.0;

  try {
    required.add_options()
        ("input,i",
         po::value<std::string>()->required(),
         "Genotype matrix file path.")
        ("ped,p",
         po::value<std::string>()->required(),
         "Path to the .ped file containing the sample phenotypes.")
        ("output,o",
         po::value<std::string>()->required(),
         "Path to output directory. Two files will be output: a simple transcript level results file, and a detailed variant level result file.");
    optional.add_options()
        ("covariates,c",
         po::value<std::string>(),
         "The covariate matrix file, tab or space separated.\nFormat = sample_id cov1 ...")
        ("bed_filter,b",
         po::value(&bed),
         "A bed file, or a comma separated list of bed files, to be used as a filter. All specified variants will be excluded.")
        ("filter,f",
         po::value<std::string>(),
         "A csv whitelist of TYPE and FUNCTION annotations. Default whitelist can be found in the filter directory.")
        ("weights,w",
         po::value(&weights),
         "A file providing weights. Replaces the CASM scores provided in the matrix file.")
        ("no_weights",
         po::bool_switch(&no_weights),
         "Disable weights.")
        ("impute_to_mean",
         po::bool_switch(&impute_to_mean),
         "Impute the to mean AF of cases for case samples, and the mean AF of controls for control samples.")
        ("whole_gene",
         po::bool_switch(&whole_gene),
         "Analyze the union of all transcripts for a gene.")
        ("nthreads,t",
         po::value<size_t>()->default_value(std::thread::hardware_concurrency() / 2 + 1),
         "The number of threads. Minimum number of threads = 2. n + 1 threads, with one parent thread and n threads processing genes.")
        ("method,m",
         po::value<std::string>()->default_value("VAAST"),
         "The statistical method to be used.\n"
         "Options: {BURDEN, CALPHA, CMC, CMC1df, RVT1, RVT2, SKAT, SKATO, SKATC, VAAST, VT, WSS}.")
        ("optimizer",
         po::value<std::string>()->default_value("irls"),
         "The optimizer used to fit the GLM.\n"
         "Options: {irls, irls_svdnewton, irls_qr, irls_qr_R, gradient_descent}.")
        ("range",
         po::value(&gene_range)->multitoken(),
         "A range of genes to analyze from the matrix file. "
         "Takes two values, a start gene number, and end gene number.\n"
         "The program will only provide results for the values in that range. "
         "Gene count starts at 1.\n"
         "Useful for starting multiple jobs on a cluster each processing part of a file.\n"
         "Note: Somewhat slower than splitting the input matrix file.")
        ("nperm",
         po::value<arma::uword>()->default_value(10000),
         "The maximum number of permutations to be performed.")
        ("mac",
         po::value<arma::uword>(),
         "Alternative or minor allele count cutoff per variant.")
        ("maf,r",
         po::value<double>()->default_value(0.5),
         "Alternative or minor allele frequency cutoff per variant. We recommend using an external sample and filtering variants based on the frequency in that sample, rather than filtering within. Can result in a reduction in power for variants near the threshold.")
        ("ma_count",
         po::bool_switch(&ma_count),
         "Change genotype matrix to minor allele counting.")
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
        ("no_detail",
         po::bool_switch(&no_detail),
         "Don't produce detailed, variant level output.")
        ("output_stats",
         po::bool_switch(&stats),
         "Write permuted statistics to .simple file following default output.")
        ("permute_out",
         po::value(&permute_set),
         "Output permutations to the given file. Exits after generating permutations.")
        ("min_minor_allele_count",
         po::value<arma::uword>()->default_value(1),
         "Minimum number of minor allele copies to test a gene.")
        ("min_variant_count",
         po::value<arma::uword>()->default_value(1),
         "Minimum number of variants to test a gene.")
        ("max_levels",
         po::value<arma::uword>()->default_value(100),
         "Maximum number of levels for a single variable. Will be split into n-1 dummy variables.")
        ("bin_epsilon",
         po::value<double>()->default_value(0.0001),
         "Odds closer together than the given value will be collapsed into a single bin for permutation.")
        ("max_perms",
         po::value(&max_perms),
         "Maximum number of permutations, used in combination with --nperm to manage memory usage. Run permutation in blocks of size nperm, up to the maximum set here. Only genes requiring additional permutation will be permuted. If you are running a small number of permutations, do not set this option.")
        ("seed",
         po::value(&seed),
         "A defined seed passed to the random number generators used for each gene.")
        ("check_testability",
         po::value(&testable),
         "Return results for genes with a minimum achievable p-value less than or equal to what is given.");
    vaast.add_options()
        ("group_size,g",
         po::value<arma::uword>()->default_value(4),
         "Group size, minor allele count threshold for grouping a variant. VAAST can collapse variants into groups of variants, dependent upon the collapse having a higher total VAAST score.")
        ("soft_maf_filter",
         po::value(&soft_maf_filter),
         "Caps the highest allele frequency for the control set in the likelihood calculation. Penalizes common variants without removing them.")
        ("biallelic",
         po::bool_switch(&biallelic),
         "Additional term for biallelic variants. For detecting potentially recessive variants.")
        ("site_penalty",
         po::value<double>()->default_value(2.0),
         "VAAST site penalty. AIC penalty applied to each site in VAAST.")
        ("alternate_grouping",
         po::bool_switch(&alternate_grouping),
         "If enabled variants are grouped all together, otherwise by VAAST 2.0 type annotation.");
    skat.add_options()
        ("kernel,k",
         po::value<std::string>()->default_value("wLinear"),
         "Kernel for use with SKAT / SKATO.\nOne of: {Linear, wLinear}.")
        ("qtl",
         po::bool_switch(&linear),
         "Analyze a quantitative trait. Values are assumed to be finite floating point values.")
        ("beta_weights",
         po::value<std::string>()->default_value("1,25"),
         "Parameters for the beta distribution. Two values, comma separated corresponding to a,b.")
        ("saddlepoint",
        po::bool_switch(&saddlepoint),
        "Force the saddlepoint approximation. Useful for highly skewed case/control sample sizes.");
    cmc.add_options()
        ("cmcmaf",
         po::value<double>()->default_value(0.005),
         "Minor allele frequency cutoff for CMC collapsing.")
        ("hotellings",
        po::bool_switch(&hotellings),
        "Use Hotellings T2 instead of a chi-square test.");
    all.add_options()
        ("wald",
         po::bool_switch(&wald),
         "Use a Wald test instead of the deviance for RVT2.")
        ("var_collapsing",
         po::bool_switch(&var_collapsing),
         "Collapse variants with <10 minor allele count into a single pseudo variant. Will convert to minor allele counting instead of alternate allele counting.")
        ("external", po::value<std::string>(), "Use external permutations.")
        ("help,h", "Print this help message.")
        ("quiet,q", "Don't print status messages.");
    hidden.add_options()
        ("nocovadj",
         po::bool_switch(&nocovadj),
         "Do Fisher-Yates shuffling instead of covariate adjusted permutation.");
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
  } catch (po::required_option &e) {
    std::cerr << "Missing required option:\n" << e.what() << "\n";
    std::cerr << visible << "\n";
    return 1;
  } catch (std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }

  const std::set<std::string> method_choices = {
      "BURDEN",
      "CALPHA",
      "CMC",
      "CMC1df",
      "VT",
      "WSS",
      "RVT1",
      "RVT2",
      "SKAT",
      "SKATO",
      "SKATC",
      "VAAST"
  };

  const std::set<std::string> kernel_choices = {
      "Linear",
      "wLinear"
  };

  const std::set<std::string> optimizer_choices = {
  	"irls",
  	"irls_svdnewton",
  	"irls_qr",
  	"irls_qr_R",
  	"gradient_descent"
  };

  if (method_choices.count(vm["method"].as<std::string>()) == 0) {
    // Method not among choices
    std::cerr
        << "Method must be one of {BURDEN, CALPHA, CMC, CMC1df, RVT1, RVT2, SKAT, SKATO, SKATC, VAAST, VT, WSS}.\n";
    std::cerr << visible << "\n";
    return 1;
  }

  if (kernel_choices.count(vm["kernel"].as<std::string>()) == 0) {
    // Method not among choices
    std::cerr << "Kernel must be one of {Linear, wLinear}.\n";
    std::cerr << visible << "\n";
    return 1;
  }

  if (optimizer_choices.count(vm["optimizer"].as<std::string>()) == 0) {
 	std::cerr << "Optimizer must be one of {irls, irls_svdnewton, irls_qr, irls_qr_R, gradient_descent}.\n";
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
  tp.nperm = vm["nperm"].as<arma::uword>();
  tp.max_perms = max_perms;

  // External permutations for Yao -- XMAT
  if (vm.count("external") > 0) {
    tp.external = true;
    tp.external_path = vm["external"].as<std::string>();
  } else {
    tp.external = false;
  }
  tp.output_stats = stats;

  tp.method = vm["method"].as<std::string>();
  tp.optimizer = vm["optimizer"].as<std::string>();
  // File paths and option status
  uint32_t pathbufsize = 1000;
  char pathbuf[pathbufsize];
  for(int i = 0; i < pathbufsize; i++) {
    pathbuf[i] = '\0';
  }
#ifdef __APPLE__
  int ret = _NSGetExecutablePath(pathbuf, &pathbufsize);
#else
  ssize_t len = readlink("/proc/self/exe", pathbuf, 1000);
#endif
  tp.program_path = pathbuf;
  ssize_t i = strlen(pathbuf);
  for (; i >= 0; i--) {
    if (pathbuf[i] == '/') {
      break;
    }
  }
  tp.program_directory = std::string(pathbuf).substr(0, i + 1);
  if (vm.count("filter") > 0) {
    tp.whitelist_path = vm["filter"].as<std::string>();
  } else {
    tp.whitelist_path = tp.program_directory + "../filter/filter_whitelist.csv";
  }
  std::cerr << "Program directory: " << tp.program_directory << std::endl;
  std::cerr << "Whitelist path: " << tp.whitelist_path << std::endl;
  tp.seed = seed;
  tp.input_path = vm["input"].as<std::string>();
  if (vm.count("covariates") > 0) {
    tp.covariates_path = vm["covariates"].as<std::string>();
  } else {
    tp.covariates_path = "";
  }
  tp.ped_path = vm["ped"].as<std::string>();
  tp.output_path = vm["output"].as<std::string>();
  tp.maf = vm["maf"].as<double>();
  tp.aaf_filter = !ma_count;
  tp.cmcmaf = vm["cmcmaf"].as<double>();
  tp.group_size = vm["group_size"].as<arma::uword>();
  tp.vaast_site_penalty = vm["site_penalty"].as<double>();
  tp.alternate_grouping = alternate_grouping;
  tp.whole_gene = whole_gene;
  tp.bed = bed;
  tp.weight = weights;
  tp.permute_set = permute_set;
  tp.bin_epsilon = vm["bin_epsilon"].as<double>();
  tp.max_levels = vm["max_levels"].as<arma::uword>();
  // Threading
  tp.nthreads = vm["nthreads"].as<size_t>();
  if (tp.nthreads < 2) {
    std::cerr << "Thread count must be >= 2." << std::endl;
    std::exit(1);
  }
  // Options
  tp.verbose = verbose;
  tp.gene_list = gene_list;
  tp.no_detail = no_detail;
  tp.top_only = top_only;
  tp.no_weights = no_weights;
  tp.impute_to_mean = impute_to_mean;
  tp.mac = vm.count("mac") > 0 ? vm["mac"].as<arma::uword>() : std::numeric_limits<unsigned long long>::max();
  tp.min_minor_allele_count = vm["min_minor_allele_count"].as<arma::uword>();
  tp.min_variant_count = vm["min_variant_count"].as<arma::uword>();
  tp.pthresh = pthresh;
  tp.soft_maf_filter = soft_maf_filter;
  // SKAT Options
  tp.kernel = vm["kernel"].as<std::string>();
  tp.qtl = linear;
  tp.saddlepoint = saddlepoint;
  tp.var_collapsing = var_collapsing;
  // Beta weights
    tp.a = std::stoi(beta_split.str(0));
    tp.b = std::stoi(beta_split.str(1));
  // Testability
  tp.testable = testable;
  tp.biallelic = biallelic;
  tp.wald = wald;

  tp.nocovadj = nocovadj || tp.covariates_path.empty();

  tp.power = false;

  // tp.alternate_permutation = tp.method == "SKATO" || tp.method == "SKAT" || tp.method == "BURDEN";
  tp.alternate_permutation = tp.nocovadj || tp.covariates_path.empty();
  tp.quantitative =
      tp.method == "RVT1" || tp.method == "RVT2" || tp.method == "SKATO" || tp.method == "SKAT" || tp.method == "BURDEN"
          || tp.method == "VT" || tp.method == "SKATC";
  tp.analytic = tp.method == "SKATO" || (tp.method == "SKAT" && tp.nperm == 0) || tp.method == "RVT1"
     || tp.method == "RVT2" || (tp.method == "CMC" && tp.nperm == 0) || (tp.method == "CMC1df" && tp.nperm == 0)
                || (tp.method == "BURDEN" && tp.nperm == 0) || tp.method == "SKATC";
  if (tp.qtl && !tp.quantitative) {
    std::cerr << "Quantitative trait analysis is only supported for the RVT1, RVT2, SKAT, SKATO, SKATC, and BURDEN methods."
              << std::endl;
    std::exit(1);
  }

  std::vector<int> range_opt;
  if (!vm["range"].empty() && (range_opt = vm["range"].as<std::vector<int>>()).size() == 2) {
    tp.range_start = range_opt[0];
    tp.range_end = range_opt[1];
  }

  if (tp.mac <= 0) {
    std::cerr << "Minor allele count cutoff must be greater than zero." << std::endl;
    std::exit(1);
  } else if (tp.mac > 500 && tp.mac < std::numeric_limits<unsigned long long>::max()) {
    std::cerr
        << "WARNING: This software is concerned with evaluating rare events. With a minor allele cutoff > 500, you should consider analyzing those variants using single marker tests."
        << std::endl;
  }

  if (tp.verbose) {
    std::cerr << "genotypes: " << tp.input_path << "\n";
    std::cerr << "covariates: " << tp.covariates_path << "\n";
    std::cerr << "ped: " << tp.ped_path << "\n";
    if (tp.bed)
      std::cerr << "bed_file: " << *tp.bed << "\n";
    if (tp.weight)
      std::cerr << "weight_file: " << *tp.weight << "\n";
    std::cerr << "output: " << tp.output_path << "\n";
    std::cerr << "method: " << tp.method << "\n";
    std::cerr << "success threshold: " << tp.success_threshold << "\n";
    std::cerr << "nthreads: " << tp.nthreads << "\n";
    std::cerr << "permutation block: " << tp.nperm << "\n";
    if (tp.max_perms)
      std::cerr << "total permutations: " << *tp.max_perms << "\n";
    if (tp.testable)
      std::cerr << "check_testability filter: " << *tp.testable << "\n";
    if (tp.maf < 0.5)
      std::cerr << "maf filter: " << tp.maf << "\n";
    if (tp.mac < std::numeric_limits<unsigned long long>::max())
      std::cerr << "mac filter: " << tp.mac << "\n";
    if (tp.gene_list)
      std::cerr << "gene list: " << *tp.gene_list << std::endl;
  }

  // Check for correct file paths
  if (!check_file_exists(tp.input_path)) {
    std::cerr << "Incorrect file path for genotypes." << std::endl;
    std::cerr << visible << "\n";
    std::exit(1);
  }
  if (vm.count("covariates") != 0 && !check_file_exists(tp.covariates_path)) {
    std::cerr << "Incorrect file path for covariates." << std::endl;
    std::cerr << visible << "\n";
    std::exit(1);
  }
  if (!check_file_exists(tp.ped_path)) {
    std::cerr << "Incorrect file path for ped." << std::endl;
    std::cerr << visible << "\n";
    std::exit(1);
  }
  if (tp.bed) {
    RJBUtil::Splitter<std::string> bed_paths(*tp.bed, ",");
    for (const auto &f : bed_paths) {
      if (!check_file_exists(std::string(f))) {
        std::cerr << "Incorrect file path for bed_file." << std::endl;
        std::cerr << visible << "\n";
        std::exit(1);
      }
    }
  }
  if (tp.weight && !check_file_exists(*tp.weight)) {
    std::cerr << "Incorrect file path for weight_file." << std::endl;
    std::cerr << visible << "\n";
    std::exit(1);
  }
  if (!check_directory_exists(tp.output_path)) {
    if (!make_directory(tp.output_path)) {
      std::cerr << "Output path is invalid. Unable to construct output directory." << std::endl;
      std::cerr << visible << "\n";
      std::exit(1);
    }
  }
  // Initialize randomization
  if(!tp.seed) {
	arma::arma_rng::set_seed_random();
  } else {
    arma::arma_rng::set_seed(*tp.seed);
  }

  if(tp.external && tp.permute_set) {
	throw std::runtime_error("Conflicting options. External permutation set will be deleted when permute_set is also included.");
  }

  if(tp.max_perms && *tp.max_perms <= 0) {
    throw std::runtime_error("If max_perms is set, it must be set to a value greater than zero.");
  }

  std::shared_ptr<Reporter> reporter = nullptr;
  reporter = std::make_shared<Reporter>(tp);
  JobDispatcher<CAPEROp, CAPERTask, Reporter> jd(tp, reporter);

  double n = timer.toc();
  std::cerr << "Elapsed time: " << n << std::endl;
  return 0;
}