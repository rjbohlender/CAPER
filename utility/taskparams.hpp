//
// Created by Bohlender,Ryan James on 10/12/18.
//

#ifndef PERMUTE_ASSOCIATE_TASKPARAMS_HPP
#define PERMUTE_ASSOCIATE_TASKPARAMS_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <string>

#include <boost/optional.hpp>

struct TaskParams {
  // Base command
  std::string base;
  std::string program_directory;
  // Permutation paramters
  arma::uword success_threshold;

  boost::optional<int> seed;

  arma::uword nperm;
  arma::uword max_levels;
  boost::optional<arma::uword> max_perms;

  bool whole_gene;

  // For external permutations
  bool external;
  std::string external_path;
  bool output_stats;

  // For SKATO, SKAT, BURDEN
  bool alternate_permutation;
  bool analytic;
  bool qtl;
  bool quantitative;
  bool saddlepoint;

  // Method
  std::string method;

  // General options
  std::string program_path;
  std::string input_path;
  std::string covariates_path;
  std::string ped_path;
  std::string whitelist_path;
  arma::uword mac;
  double maf;
  double min_variant_count;
  double min_minor_allele_count;
  bool nocovadj;
  bool no_weights;
  bool impute_to_mean;
  bool aaf_filter;

  boost::optional<int> range_start;
  boost::optional<int> range_end;

  boost::optional<std::string> bed;
  boost::optional<std::string> weight;

  double bin_epsilon;

  double soft_maf_filter;
  double vaast_site_penalty;
  bool alternate_grouping;

  bool verbose;

  // Output permutations
  boost::optional<std::string> permute_set;

  // Alternate p-value threshold
  boost::optional<double> pthresh;

  // Choose optimizer implementation for GLM
  std::string optimizer;

  // Detailed VAAST output
  std::string full_command;
  std::string output_path;
  bool top_only;

  // VAAST
  arma::uword group_size;
  boost::optional<double> testable;
  bool biallelic;

  // CMC
  double cmcmaf;

  size_t nthreads;

  // Gene list
  boost::optional<std::string> gene_list;
  bool no_detail;

  // Run power analysis
  bool power;
  arma::vec alpha;
  arma::uword bootstrap_reps;
  std::vector<arma::uword> ncases;
  std::vector<arma::uword> ncontrols;

  // SKAT Parameters
  std::string kernel; // Kernel selection
  int a; // Beta weight parameters
  int b; // Beta weight parameters
};

#endif //PERMUTE_ASSOCIATE_TASKPARAMS_HPP
