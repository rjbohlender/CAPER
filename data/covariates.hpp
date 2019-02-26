//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_COVARIATES_HPP
#define PERMUTE_ASSOCIATE_COVARIATES_HPP

#define ARMA_DONT_USE_WRAPPER

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <armadillo>

#include "../utility/split.hpp"
#include "../statistics/logistic_regression.hpp"
#include "permutation.hpp"

class Covariates {
public:
  explicit Covariates(const std::string& ifile, const std::string& pedfile, bool linear=false);
  explicit Covariates(std::stringstream& ss);

  Covariates(const Covariates &cov) = default;

  Covariates(Covariates &&cov) = delete; // Covariates should never be moved and always copied.

  Covariates &operator=(const Covariates &cov) = default;

  void print();

  // Getters and Setters
  arma::colvec &get_phenotype_vector();
  void set_phenotype_vector(arma::colvec &vec);
  void set_phenotype_vector(std::vector<int32_t> &vec);

  arma::uword get_nsamples() const;
  arma::uword get_ncases() const;

  arma::mat &get_covariate_matrix();
  arma::colvec &get_odds();
  arma::colvec &get_original_phenotypes();
  arma::vec &get_fitted();
  arma::uvec &get_indices();
  arma::vec get_residuals() const;
  arma::vec &get_mean();
  arma::rowvec &get_coef();

  // Permute
  void shuffle();

  // Free memory
  void clear();

  // Sort covariates
  void sort_covariates(std::string &header);
  bool is_sorted();

private:
  CRandomMersenne crand;
  std::vector<std::string> samples_;
  arma::colvec phenotypes_; // Possibly permuted phenotype vector
  arma::colvec original_; // Original phenotype vector
  arma::mat design_; // Design matrix
  arma::colvec odds_;
  arma::vec fitted_;
  arma::uvec indices_;
  arma::vec mean_; // Mean of MFNCH
  arma::vec eta_;
  arma::rowvec coef_;

  bool sorted_;
  bool linear_;

  unsigned long nsamples_;
  unsigned long ncases_;

  void parse(const std::string& ifile, const std::string& pedfile);
  void parse(std::stringstream& ss);

  void fit_null();
};

#endif //PERMUTE_ASSOCIATE_COVARIATES_HPP
