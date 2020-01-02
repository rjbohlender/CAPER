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
  void set_phenotype_vector(arma::vec vec);
  void set_phenotype_vector(std::vector<int32_t> vec);

  arma::uword get_nsamples() const;
  arma::uword get_ncases() const;

  arma::mat &get_covariate_matrix();
  arma::vec &get_odds();
  arma::vec &get_original_phenotypes();
  arma::vec &get_fitted();
  arma::vec get_residuals() const;
  arma::vec &get_mean();
  arma::rowvec &get_coef();

  // Permute
  void refit_permuted();

  // Free memory
  void clear();

  // Sort covariates
  void sort_covariates(std::string &header);
  bool is_sorted();

private:
  unsigned long nsamples_;
  unsigned long ncases_;

  CRandomMersenne crand;
  std::vector<std::string> cov_samples_;
  std::vector<std::string> ped_samples_;
  arma::vec phenotypes_; // Possibly permuted phenotype vector
  arma::vec original_; // Original phenotype vector
  arma::mat design_; // Design matrix
  arma::vec odds_;
  arma::vec fitted_;
  arma::vec mean_; // Mean of MFNCH
  arma::vec eta_;
  arma::rowvec coef_;

  // Permutation
  arma::vec p_odds_;
  arma::vec p_fitted_;
  arma::vec p_eta_;
  arma::rowvec p_coef_;

  bool sorted_;
  bool linear_;

  void parse(const std::string& ifile, const std::string& pedfile);
  void parse(std::stringstream& ss);

  void fit_null();
};

#endif //PERMUTE_ASSOCIATE_COVARIATES_HPP
