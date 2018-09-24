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

#include "../utility/split.hpp"
#include "weight.hpp"

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

  arma::vec &get_weights(const std::string &k);
  void set_weights(const std::string &k, arma::vec &weights);

  bool is_weighted(const std::string &k);

  void clear();

private:
  std::map<std::string, bool> weights_set_;

  unsigned long nsamples_;
  std::map<std::string, unsigned long> nvariants_;
  std::string gene_;
  std::vector<std::string> transcripts_;
  std::map<std::string, std::vector<std::string>> positions_;

  std::map<std::string, arma::mat> genotypes_;
  std::map<std::string, arma::vec> weights_;

  std::string header_;

  void parse(std::stringstream &ss);
};

#endif //PERMUTE_ASSOCIATE_GENE_HPP
