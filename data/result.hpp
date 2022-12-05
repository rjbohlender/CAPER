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
#include "../utility/taskparams.hpp"

struct Result {
  Result();
  Result(std::string gene, std::string transcript, bool skippable);
  Result(const Result &res) = default;
  Result(Result &&res) noexcept;

  Result &operator=(const Result &rhs) = default;
  Result &operator=(Result &&rhs) noexcept;

  friend std::ostream &operator<<(std::ostream &stream, Result &rhs);

  Result &combine(const Result &res, const TaskParams &tp);
  void update_ci();
  void calc_exact_p();
  void calc_exact_p(double n1, double n2);

  bool is_set = false;
  std::string gene;
  std::string transcript;
  std::pair<double, double> empirical_ci;
  std::pair<double, double> empirical_midci;
  size_t successes;
  double mid_successes;
  size_t permutations;
  size_t min_success_at;
  double rand_perms;
  double original;
  double empirical_p;
  double empirical_midp;
  double exact_p;
  double mgit_p;
  size_t mgit_successes;
  bool done;
  bool skippable;
  std::vector<double> permuted;
  int rank;
  bool testable;
  double odds;
  double nmac, nmaj;
  // std::map<double, double> permuted;

  // Output permuted statistics
  bool output_stats;

  // Odds Ratio specific
  double or_p;
  double case_alt, case_ref;
  double cont_alt, cont_ref;
};


#endif //PERMUTE_ASSOCIATE_RESULT_HPP
