//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_GENE_HPP
#define PERMUTE_ASSOCIATE_GENE_HPP

#define ARMA_DONT_USE_WRAPPER

#include "../utility/split.hpp"
#include "../utility/taskparams.hpp"
#include "covariates.hpp"
#include "filter.hpp"
#include "matrix_indices.hpp"
#include "result.hpp"
#include "weight.hpp"

#include <armadillo>
#include <iostream>
#include <iterator>
#include <unordered_map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <stack>

void print_comma_sep(std::vector<std::string> &x, std::ostream &os);
void print_comma_sep(const std::vector<std::string> &x, std::ostream &os);
void print_semicolon_sep(arma::uvec &x, std::ostream &os);

/**
 * @brief Container for single gene multiple transcript data
 */
class Gene {
public:
  std::string gene_name;

  std::vector<std::string> transcripts;
  std::vector<arma::uword> columns;
  std::vector<std::string> samples;

  std::unordered_map<std::string, std::vector<std::string>> positions;
  std::unordered_map<std::string, arma::sp_mat> genotypes;
  std::unordered_map<std::string, arma::vec> weights;
  std::unordered_map<std::string, std::vector<std::string>> function;
  std::unordered_map<std::string, std::vector<std::string>> annotation;
  std::unordered_map<std::string, std::vector<std::string>> reference;
  std::unordered_map<std::string, std::vector<std::string>> alternate;
  std::unordered_map<std::string, std::vector<std::string>> type;
  std::unordered_map<std::string, std::stack<int>> to_remove;

  bool testable;
  bool skippable;

  Gene(std::stringstream &ss, const std::shared_ptr<Covariates> &cov,
       unsigned long nsamples,
       std::unordered_map<std::string, arma::uword> nvariants,
       const Weight &weight, const TaskParams &tp, Filter &filter);

  void print();

  void set_matrix(const std::string &k, arma::sp_mat &data);
  void set_matrix(const std::string &k, arma::sp_mat &&data);
  std::vector<std::string> &get_transcripts();
  std::vector<std::string> &get_positions(const std::string &k);
  [[nodiscard]] const std::vector<std::string> &get_samples() const;
  arma::vec &get_weights(const std::string &k);
  void set_weights(const std::string &k, arma::vec &new_weights);
  void set_scores(const std::string &k, arma::vec &scores);
  [[nodiscard]] std::string get_detail() const;
  [[nodiscard]] std::unordered_map<std::string, std::string> get_vaast() const;
  [[nodiscard]] bool is_skippable() const;
  bool is_polymorphic(const std::string &k);
  void generate_detail(Covariates &cov,
                       std::unordered_map<std::string, Result> &results);
  void generate_vaast(Covariates &cov);
  void clear(Covariates &cov, std::unordered_map<std::string, Result> &results,
             const TaskParams &tp);

  static std::string form_variant_id(RJBUtil::Splitter<std::string> &splitter);

private:
  uint64_t nsamples; // Number of used samples, <= total_samples
  uint64_t total_samples; // Number of samples in the matrix header
  std::unordered_map<std::string, arma::uword> nvariants;

  std::unordered_map<std::string, arma::vec> variant_scores; // Stored if detail is true
  std::unordered_map<std::string, bool> polymorphic;

  const TaskParams &tp;

  std::string detail_;
  std::unordered_map<std::string, std::string> vaast_;

  std::unordered_map<std::string, arma::sp_mat> missing_variant_carriers_;

  void parse(std::stringstream &ss, const std::shared_ptr<Covariates> &cov,
             Filter &filter);
  static std::string compress_adjacent(arma::uvec &samples);

  bool check_testability(const std::string &transcript, Covariates &cov,
                const std::vector<double> &permuted);

  std::stringstream transcript_union(std::stringstream &ss,
                                     const std::shared_ptr<Covariates> &cov,
                                     Filter &filter);

  void aaf_filter();
  void maf_filter();
  void collapse_variants();
  void impute_to_mean(const std::shared_ptr<Covariates> &cov);
  void update_weights(const Weight &weight);
};

#endif // PERMUTE_ASSOCIATE_GENE_HPP
