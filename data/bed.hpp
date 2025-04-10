//
// Created by Bohlender,Ryan James on 8/9/18.
//

#ifndef PERMUTE_ASSOCIATE_BED_HPP
#define PERMUTE_ASSOCIATE_BED_HPP

#include <algorithm>
#include <exception>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>

#include "../utility/split.hpp"

/**
 * @brief Convenience wrapper class for position ranges
 */
struct BedRange {
  explicit BedRange(std::pair<int, int> &range);
  BedRange(int first, int second);

  std::pair<int, int> range;

  // Ordering
  bool operator<(int rhs) const;
  bool operator>(int rhs) const;
  bool operator==(int rhs) const;
  bool operator<=(int rhs) const;
  bool operator>=(int rhs) const;

  // Arithmetic
  int operator-(int rhs) const;
};

class Bed {
public:
  // Constructors
  Bed() = default;
  explicit Bed(const std::string &ifile);
  explicit Bed(std::stringstream &iss);

  bool check_variant(const std::string &variant);

  bool empty();
  unsigned long size();

private:
  std::map<std::string, std::vector<BedRange>> ranges_;
  std::unordered_set<std::string> variants_;
};

#endif // PERMUTE_ASSOCIATE_BED_HPP
