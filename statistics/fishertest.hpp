//
// Created by Bohlender,Ryan James on 2019-06-06.
//

#ifndef PERMUTE_ASSOCIATE_FISHERTEST_HPP
#define PERMUTE_ASSOCIATE_FISHERTEST_HPP

#include "../data/gene.hpp"
#include "../data/covariates.hpp"

class FisherTest {
public:
  double case_ref, case_alt;
  double cont_ref, cont_alt;
  // Estimate odds ratio via Fisher's Exact Test
  FisherTest(Gene &gene, arma::vec &Y, const std::string &ts);

  auto get_or() -> double; // Return the OR estimate
  auto get_pval() -> double; // Return the P value estimate

private:
  double or_;
  double p_;
};

class CAFisherTest {
public:
  CAFisherTest(Gene &gene, arma::vec &Y, const std::string &ts);

  auto get_or() -> double; // Return the OR estimate
  auto get_pval() -> double; // Return the P value estimate

private:
  CMultiFishersNCHypergeometric cmfnch_;

};
#endif //PERMUTE_ASSOCIATE_FISHERTEST_HPP
