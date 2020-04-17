//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_GENE_HPP
#define PERMUTE_ASSOCIATE_GENE_HPP

#define ARMA_DONT_USE_WRAPPER

#include <string>
#include <sstream>
#include <iostream>
#include <armadillo>
#include <iterator>
#include <map>
#include <unordered_map>

#include "../utility/split.hpp"
#include "weight.hpp"
#include "covariates.hpp"
#include "result.hpp"
#include "../utility/taskparams.hpp"

arma::uvec setdiff(arma::uvec x, arma::uvec y);

void print_comma_sep(arma::uvec &x, std::ostream &os);
void print_comma_sep(std::vector<std::string> &x, std::ostream &os);
void print_comma_sep(const std::vector<std::string> &x, std::ostream &os);
void print_semicolon_sep(arma::uvec &x, std::ostream &os);

/**
 * @brief Container for single gene multiple transcript data
 */
class Gene {
public:
  Gene(std::stringstream &ss, unsigned long nsamples, std::map<std::string, arma::uword> &nvariants,
       const Weight &weight, TaskParams tp, arma::vec &phenotypes);

  void print();

  // Member access
  std::string &get_gene();
  arma::sp_mat &get_matrix(const std::string &k);
  arma::sp_mat &get_missing(const std::string &k);
  void set_matrix(const std::string &k, arma::sp_mat &&data);
  std::vector<std::string> &get_transcripts();
  arma::uword get_nvariants(const std::string &k);
  std::vector<std::string> &get_positions(const std::string &k);
  std::vector<std::string> &get_samples();
  arma::vec &get_weights(const std::string &k);
  void set_weights(const std::string &k, arma::vec &weights);
  arma::vec &get_scores(const std::string &k);
  void set_scores(const std::string &k, arma::vec &scores);
  auto get_detail() -> std::string;
  auto get_vaast() -> std::map<std::string, std::string>;
  auto is_skippable() const -> bool;
  auto is_polymorphic(const std::string &k) -> bool;
  auto is_weighted(const std::string &k) -> bool;
  auto is_testable() const -> bool;
  auto generate_detail(Covariates &cov, std::unordered_map<std::string, Result> &results, TaskParams &tp) -> void;
  auto generate_vaast(Covariates &cov) -> void;
  auto clear(Covariates &cov, std::unordered_map<std::string, Result> &results, TaskParams &tp) -> void;

private:
  std::map<std::string, bool> weights_set_;

  unsigned long nsamples_;
  std::map<std::string, arma::uword> nvariants_;
  std::string gene_;
  std::vector<std::string> transcripts_;
  std::map<std::string, std::vector<std::string>> positions_;

  std::map<std::string, arma::sp_mat> genotypes_;
  std::map<std::string, arma::vec> weights_;
  std::vector<std::string> samples_;

  std::map<std::string, arma::vec> variant_scores_; // Stored if detail is true
  std::map<std::string, double> odds_;
  bool testable_;
  bool skippable_;
  std::map<std::string, bool> polymorphic_;

  TaskParams tp_;

  std::string header_;
  std::string detail_;
  std::map<std::string, std::string> vaast_;

  std::map<std::string, arma::sp_mat> missing_variant_carriers_;

  void parse(std::stringstream &ss, arma::vec &phenotypes);

  auto testable(const std::string &k, Covariates &cov, TaskParams &tp) -> bool;
};

#endif //PERMUTE_ASSOCIATE_GENE_HPP
