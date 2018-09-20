//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_COVARIATES_HPP
#define PERMUTE_ASSOCIATE_COVARIATES_HPP

#include "../utility/split.hpp"
#include "../statistics/logistic_regression.hpp"
#include "permutation.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <armadillo>

class Covariates {
public:
  explicit Covariates(const std::string& ifile, const std::string& pedfile);
  explicit Covariates(std::stringstream& ss);

  Covariates(const Covariates &cov) = default;

  Covariates(Covariates &&cov) = delete; // Covariates should never be moved and always copied.

  Covariates &operator=(const Covariates &cov) = default;

  void print();

  // Getters and Setters
  arma::colvec &get_phenotype_vector();
  void set_phenotype_vector(arma::colvec &vec);
  void set_phenotype_vector(std::vector<int32_t> &vec);

  unsigned long get_nsamples();
  unsigned long get_ncases();

  arma::mat &get_covariate_matrix();
  arma::colvec &get_odds();
  arma::colvec &get_original_phenotypes();
  arma::vec &get_probability();
  arma::uvec &get_indices();
  arma::vec &get_mean();

  // Permute
  void shuffle();

  // Free memory
  void clear();

  // Sort covariates
  void sort_covariates(std::string &header);
  bool is_sorted();

private:
  std::vector<std::string> samples_;
  arma::colvec phenotypes_;
  arma::colvec original_;
  arma::mat covariates_;
  arma::colvec odds_;
  arma::vec prob_;
  arma::uvec indices_;
  arma::vec mean_; // Mean of MFNCH
  arma::vec eta_;

  bool sorted_;

  unsigned long nsamples_;
  unsigned long ncases_;

  void parse(const std::string& ifile, const std::string& pedfile);
  void parse(std::stringstream& ss);

  void calculate_odds();
};

#endif //PERMUTE_ASSOCIATE_COVARIATES_HPP
