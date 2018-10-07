//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_GENE_HPP
#define PERMUTE_ASSOCIATE_GENE_HPP

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

arma::uvec setdiff(arma::uvec x, arma::uvec y);
void print_comma_sep(arma::uvec &x, std::ostream &os);
void print_comma_sep(std::vector<std::string> &x, std::ostream &os);
void print_comma_sep(const std::vector<std::string> &x, std::ostream &os);
void print_semicolon_sep(arma::uvec &x, std::ostream &os);

class Gene {
public:
  Gene(std::stringstream &ss, unsigned long nsamples, std::map<std::string, unsigned long> &nvariants, const Weight &weight);

  void print();

  // Getters
  std::string &get_gene();
  arma::mat &get_matrix(const std::string &k);
  std::vector<std::string> &get_transcripts();
  unsigned long get_nvariants(const std::string &k);
  std::vector<std::string> &get_positions(const std::string &k);
  std::vector<std::string> &get_samples();

  arma::vec &get_weights(const std::string &k);
  void set_weights(const std::string &k, arma::vec &weights);

  arma::vec &get_scores(const std::string &k);
  void set_scores(const std::string &k, arma::vec &scores);

  std::string get_detail();

  bool is_weighted(const std::string &k);

  void generate_detail(Covariates &cov, std::unordered_map<std::string, Result> &results);

  void clear(Covariates &cov, std::unordered_map<std::string, Result> &results);

private:
  std::map<std::string, bool> weights_set_;

  unsigned long nsamples_;
  std::map<std::string, unsigned long> nvariants_;
  std::string gene_;
  std::vector<std::string> transcripts_;
  std::map<std::string, std::vector<std::string>> positions_;

  std::map<std::string, arma::mat> genotypes_;
  std::map<std::string, arma::vec> weights_;
  std::vector<std::string> samples_;

  std::map<std::string, arma::vec> variant_scores_; // Stored if detail is true

  std::string header_;
  std::string detail_;

  void parse(std::stringstream &ss);
};

#endif //PERMUTE_ASSOCIATE_GENE_HPP
