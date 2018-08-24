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
#include "casm.hpp"

class Gene {
public:
  Gene(std::stringstream &ss, unsigned long nsamples, std::map<std::string, unsigned long> &nvariants, const CASM &casm);

  void print();

  // Getters
  arma::mat &get_matrix(const std::string &k);
  std::string &get_gene();
  std::vector<std::string> &get_transcripts();
  arma::vec &get_weights(const std::string &k);
  unsigned long get_nvariants(const std::string &k);
  std::vector<std::string> &get_positions(const std::string &k);

  void set_weights(const std::string &k, arma::vec &weights);

  void clear();

private:
  unsigned long nsamples_;
  std::map<std::string, unsigned long> nvariants_;
  std::string gene_;
  std::vector<std::string> transcripts_;
  std::map<std::string, std::vector<std::string>> positions_;

  std::map<std::string, arma::mat> genotypes_;
  std::map<std::string, arma::vec> weights_;

  // TODO Strip out header
  std::string header_;

  void parse(std::stringstream &ss);
};

#endif //PERMUTE_ASSOCIATE_GENE_HPP
