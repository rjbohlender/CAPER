//
// Created by Bohlender,Ryan James on 8/14/18.
//

#ifndef PERMUTE_ASSOCIATE_RESULT_HPP
#define PERMUTE_ASSOCIATE_RESULT_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <vector>


struct Result {
  Result();
  Result(const std::string &gene, const std::string &transcript);
  Result(const Result &res) = default;
  Result(Result &&res) noexcept;

  Result &operator=(const Result &rhs) = default;
  Result &operator=(Result &&rhs) noexcept;

  friend std::ostream &operator<<(std::ostream &stream, const Result &rhs);

  Result &combine(const Result& res);
  void set_rank(int rank);

  std::string gene;
  std::string transcript;
  int successes;
  double mid_successes;
  int permutations;
  int min_success_at;
  double original;
  double empirical_p;
  double empirical_midp;
  double mgit_p;
  int mgit_successes;
  bool done;
  std::vector<double> permuted;
  int rank;
  // std::map<double, double> permuted;
};


#endif //PERMUTE_ASSOCIATE_RESULT_HPP
