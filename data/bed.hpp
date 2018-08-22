//
// Created by Bohlender,Ryan James on 8/9/18.
//

#ifndef PERMUTE_ASSOCIATE_BED_HPP
#define PERMUTE_ASSOCIATE_BED_HPP

#include <string>
#include <map>
#include <utility>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <exception>

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
  explicit Bed(std::stringstream &ss);

  bool check_variant(const std::string &chr, const std::string &pos);
  bool check_variant(const std::string &chr, int pos);

  bool empty();
  unsigned long size();

  // Number of entries in ranges_ for a given chromosome
  unsigned long chromosome_count(const std::string &k);
private:
  std::map<std::string, std::vector<BedRange>> ranges_;
};

#endif //PERMUTE_ASSOCIATE_BED_HPP
