//
// Created by Bohlender,Ryan James on 10/12/18.
//

#ifndef PERMUTE_ASSOCIATE_TASKPARAMS_HPP
#define PERMUTE_ASSOCIATE_TASKPARAMS_HPP

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif

#include <armadillo>
#include <string>

#include <boost/optional.hpp>

struct TaskParams {
  // Base command
  std::string base;
  std::string program_directory;
  // Permutation paramters
  arma::uword success_threshold = 0;

  boost::optional<int> seed;

  arma::uword nperm = 0;
  arma::uword max_levels = 100;
  boost::optional<arma::uword> max_perms;

  bool whole_gene = false;

  // For external permutations
  bool external = false;
  std::string external_path;
  bool output_stats = false;

  // For SKATO, SKAT, BURDEN
  bool alternate_permutation = false;
  bool analytic = false;
  bool qtl = false;
  bool quantitative = false;
  bool saddlepoint = false;

  // Method
  std::string method;

  // General options
  std::string program_path;
  std::string input_path;
  std::string covariates_path;
  std::string ped_path;
  std::string whitelist_path;
  arma::uword mac = 0;
  double maf = 0.5;
  double min_variant_count = 1;
  double min_minor_allele_count = 1;
  bool nocovadj = false;
  bool no_weights = false;
  bool impute_to_mean = false;
  bool aaf_filter = false;

  boost::optional<int> range_start;
  boost::optional<int> range_end;

  boost::optional<std::string> bed;
  boost::optional<std::string> weight;

  double bin_epsilon = 0.0001;

  double soft_maf_filter = 0.0;
  double vaast_site_penalty = 0.0;
  bool alternate_grouping = false;

  bool verbose = false;

  // RVT test options
  bool wald = false;

  // Output permutations
  boost::optional<std::string> permute_set;

  // Alternate p-value threshold
  boost::optional<double> pthresh;

  // Choose optimizer implementation for GLM
  std::string optimizer;

  // Detailed VAAST output
  std::string full_command;
  std::string output_path;
  bool top_only = false;

  // VAAST
  arma::uword group_size = 4;
  boost::optional<double> testable;
  bool biallelic = false;

  // CMC
  double cmcmaf = 0.005;
  bool hotellings = false;

  size_t nthreads = 1;

  // Gene list
  boost::optional<std::string> gene_list;
  bool no_detail = false;

  // Run power analysis
  bool power = false;
  arma::vec alpha;
  arma::uword bootstrap_reps = 0;
  std::vector<arma::uword> ncases;
  std::vector<arma::uword> ncontrols;

  // SKAT Parameters
  std::string kernel; // Kernel selection
  int a = 1; // Beta weight parameters
  int b = 25; // Beta weight parameters
  bool var_collapsing = false;
};

#endif //PERMUTE_ASSOCIATE_TASKPARAMS_HPP
