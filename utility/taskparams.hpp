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
  // Permutation paramters
  arma::uword success_threshold;

  arma::uword stage_1_permutations;
  arma::uword stage_2_permutations;
  arma::uword total_permutations;

  // For external permutations
  bool external;
  std::string external_path;
  bool output_stats;

  // For SKATO, SKAT, BURDEN
  bool alternate_permutation;
  bool analytic;
  bool linear;
  bool quantitative;

  // Method
  std::string method;

  // General options
  std::string program_path;
  std::string genotypes_path;
  std::string covariates_path;
  std::string ped_path;
  arma::uword mac;
  double maf;
  bool covadj;
  bool cov_adjusted; // Methods that cannot be covariate adjusted permute a subset

  boost::optional<arma::uword> approximate;
  arma::uword maj_nbins;
  boost::optional<int> range_start;
  boost::optional<int> range_end;

  boost::optional<std::string> bed;
  boost::optional<std::string> weight;

  bool verbose;

  // Output permutations
  boost::optional<std::string> permute_set;

  // Alternate p-value threshold
  boost::optional<double> pthresh;

  // Detailed VAAST output
  std::string full_command;
  std::string output_path;
  bool top_only;

  // VAAST
  arma::uword group_size;
  bool score_only_minor;
  bool score_only_alternative;
  bool testable;
  bool biallelic;

  // CMC
  double cmcmaf;

  size_t nthreads;

  // Gene list
  boost::optional<std::string> gene_list;
  bool nodetail;

  // Run power analysis
  bool power;
  arma::vec alpha;
  arma::uword bootstrap_reps;
  std::vector<arma::uword> ncases;
  std::vector<arma::uword> ncontrols;

  // SKAT Parameters
  std::string kernel; // Kernel selection
  bool adjust; // Sample size adjustment
  int a; // Beta weight parameters
  int b; // Beta weight parameters
};

#endif //PERMUTE_ASSOCIATE_TASKPARAMS_HPP
