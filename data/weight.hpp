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

class Weight {
public:
  Weight() = default;
  explicit Weight(const std::string &ifile);
  explicit Weight(std::stringstream &ifile);

  double get(const std::string &k) const;
  double get(const std::string &k);

  bool empty() const;
  bool empty();

private:
  std::unordered_map<std::string, double> scores_;

};

#endif //PERMUTE_ASSOCIATE_CASM_HPP
