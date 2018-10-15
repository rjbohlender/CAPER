//
// Created by Bohlender,Ryan James on 10/12/18.
//

#ifndef PERMUTE_ASSOCIATE_TASKPARAMS_HPP
#define PERMUTE_ASSOCIATE_TASKPARAMS_HPP

#include <armadillo>
#include <string>

#include <boost/optional.hpp>

struct TaskParams {
  // Permutation paramters
  arma::uword success_threshold;

  arma::uword stage_1_permutations;
  arma::uword stage_2_permutations;
  arma::uword total_permutations;

  // For SKATO, SKAT, BURDEN
  bool alternate_permutation;

  // Method
  std::string method;

  // General options
  std::string program_path;
  std::string genotypes_path;
  std::string covariates_path;
  std::string ped_path;

  boost::optional<std::string> bed;
  boost::optional<std::string> weight;

  bool verbose;

  // Output permutations
  boost::optional<std::string> permute_set;

  // Detailed VAAST output
  std::string full_command;
  std::string output_path;

  // VAAST
  arma::uword group_size;
  bool score_only_minor;
  bool score_only_alternative;
  bool testable;

  // CMC
  double maf;

  size_t nthreads;

  // Gene list
  boost::optional<std::string> gene_list;

  // SKAT Parameters
  std::string kernel; // Kernel selection
  bool adjust; // Sample size adjustment
  double a; // Beta weight parameters
  double b; // Beta weight parameters
};

#endif //PERMUTE_ASSOCIATE_TASKPARAMS_HPP
