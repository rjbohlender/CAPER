//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_GENE_HPP
#define PERMUTE_ASSOCIATE_GENE_HPP

#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>

#include "../utility/split.hpp"
#include "../utility/taskparams.hpp"
#include "covariates.hpp"
#include "filter.hpp"
#include "matrix_indices.hpp"
#include "result.hpp"
#include "weight.hpp"

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
  std::string gene_name;

  Gene(std::stringstream &ss, std::shared_ptr<Covariates> cov,
       unsigned long nsamples, std::map<std::string, arma::uword> &nvariants,
       const Weight &weight, TaskParams tp, Filter &filter);

  void print();

  arma::sp_mat &get_matrix(const std::string &k);
  arma::sp_mat &get_missing(const std::string &k);
  void set_matrix(const std::string &k, arma::sp_mat &data);
  void set_matrix(const std::string &k, arma::sp_mat &&data);
  std::vector<std::string> &get_transcripts();
  arma::uword get_nvariants(const std::string &k);
  std::vector<std::string> &get_positions(const std::string &k);
  const std::vector<std::string> &get_samples() const;
  arma::vec &get_weights(const std::string &k);
  void set_weights(const std::string &k, arma::vec &weights);
  arma::vec &get_scores(const std::string &k);
  void set_scores(const std::string &k, arma::vec &scores);
  std::string get_detail() const;
  std::map<std::string, std::string> get_vaast() const;
  bool is_skippable() const;
  bool is_polymorphic(const std::string &k);
  bool is_testable() const;
  void generate_detail(Covariates &cov,
                       std::unordered_map<std::string, Result> &results,
                       const TaskParams &tp);
  void generate_vaast(Covariates &cov);
  void clear(Covariates &cov, std::unordered_map<std::string, Result> &results,
             const TaskParams &tp);

  static std::string form_variant_id(RJBUtil::Splitter<std::string> &splitter);

private:
  unsigned long nsamples_;
  std::map<std::string, arma::uword> nvariants_;
  std::vector<std::string> transcripts_;
  std::map<std::string, std::vector<std::string>> positions_;
  std::vector<arma::uword> columns_;

  std::map<std::string, arma::sp_mat> genotypes_;
  std::map<std::string, arma::vec> weights_;
  std::map<std::string, std::vector<std::string>> function_;
  std::map<std::string, std::vector<std::string>> annotation_;
  std::map<std::string, std::vector<std::string>> reference_;
  std::map<std::string, std::vector<std::string>> alternate_;
  std::map<std::string, std::vector<std::string>> type_;
  std::vector<std::string> samples_;

  std::map<std::string, arma::vec> variant_scores_; // Stored if detail is true
  bool testable_;
  bool skippable_;
  std::map<std::string, bool> polymorphic_;

  TaskParams tp_;

  std::string detail_;
  std::map<std::string, std::string> vaast_;

  std::map<std::string, arma::sp_mat> missing_variant_carriers_;

  void parse(std::stringstream &ss, const std::shared_ptr<Covariates> &cov,
             Filter &filter);
  std::string compress_adjacent(arma::uvec &samples);

  bool testable(const std::string &transcript, Covariates &cov,
                const std::vector<double> &permuted);

  std::stringstream transcript_union(std::stringstream &ss,
                                     const std::shared_ptr<Covariates> &cov,
                                     Filter &filter);
};

#endif // PERMUTE_ASSOCIATE_GENE_HPP
