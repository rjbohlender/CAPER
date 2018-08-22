//
// Created by Bohlender,Ryan James on 8/14/18.
//

#ifndef PERMUTE_ASSOCIATE_CASM_HPP
#define PERMUTE_ASSOCIATE_CASM_HPP

#include <string>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>

#include "../utility/split.hpp"

class CASM {
public:
  CASM() = default;
  explicit CASM(const std::string &ifile);
  explicit CASM(std::stringstream &iss, bool log_vals); // For testing

  double get(const std::string &k) const;
  double get(const std::string &k);

  bool empty() const;
  bool empty();

private:
  std::unordered_map<std::string, double> scores_;

};

#endif //PERMUTE_ASSOCIATE_CASM_HPP
