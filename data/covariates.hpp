//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_COVARIATES_HPP
#define PERMUTE_ASSOCIATE_COVARIATES_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../utility/taskparams.hpp"
#include "permutation.hpp"

class Covariates {
public:
  explicit Covariates(TaskParams tp);
  Covariates(std::stringstream &ped_ss, std::stringstream &cov_ss,
             TaskParams tp);

  Covariates(const Covariates &cov) = default;

  // Covariates(Covariates &&cov) = delete; // Covariates should never be moved and always copied.

  Covariates &operator=(const Covariates &cov) = default;

  void print();

  // Getters and Setters
  arma::colvec &get_phenotype_vector();
  void set_phenotype_vector(const arma::vec &vec);
  void set_phenotype_vector(const std::vector<int32_t> &vec);
  bool contains(const std::string &sample) const;
  bool contains(const std::string_view &sample) const;

  arma::uword get_nsamples() const;
  arma::uword get_ncases() const;

  arma::mat &get_covariate_matrix();
  arma::vec &get_odds();
  arma::vec &get_original_phenotypes();
  arma::vec &get_fitted();
  arma::vec get_residuals() const;
  arma::vec &get_mean();
  arma::vec &get_coef();
  std::vector<std::string> get_samples();

  // Permute
  void refit_permuted();

  // Free memory
  void clear();

  // Sort covariates
  void sort_covariates(const std::string &header);
  bool is_sorted() const;

private:
  TaskParams tp_;
  unsigned long nsamples_;
  unsigned long ncases_;
  unsigned long ncontrols_;

  CRandomMersenne crand;
  std::vector<std::string> cov_samples_;
  std::vector<std::string> ped_samples_ordered_;
  std::unordered_set<std::string> ped_samples_;
  std::unordered_set<std::string> skip_;  // Skip samples with missing cov values
  std::unordered_map<std::string, double> sample_phen_map_;
  arma::vec phenotypes_; // Possibly permuted phenotype vector
  arma::vec original_; // Original phenotype vector
  arma::mat design_; // Design matrix
  arma::vec odds_;
  arma::vec fitted_;
  arma::vec mean_; // Mean of MFNCH
  arma::vec eta_;
  arma::vec coef_;

  // Permutation
  arma::vec p_odds_;
  arma::vec p_fitted_;
  arma::vec p_eta_;
  arma::vec p_coef_;

  bool sorted_;
  bool linear_;

  void parse_ped(const std::string& pedfile, bool cov_provided);
  void parse_cov(const std::string &covfile);
  void parse(const std::string& covfile, const std::string& pedfile);
  void parse(std::stringstream &ped_ss, std::stringstream &cov_ss);

  void fit_null();

  static void sort_by_index(std::vector<std::string> &samples, const arma::uvec &indices);
};

#endif //PERMUTE_ASSOCIATE_COVARIATES_HPP
