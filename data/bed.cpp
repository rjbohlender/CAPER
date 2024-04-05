//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include <iostream>

#include "../utility/filesystem.hpp"
#include "../utility/filevalidator.hpp"
#include "bed.hpp"

BedRange::BedRange(std::pair<int, int> &range) : range(range) {}

BedRange::BedRange(int first, int second) : range({first, second}) {}

bool BedRange::operator<(int rhs) const { return range.second < rhs; }
bool BedRange::operator>(int rhs) const { return range.first > rhs; }
bool BedRange::operator==(int rhs) const {
  return range.first <= rhs && range.second >= rhs;
}
bool BedRange::operator<=(int rhs) const { return range.first <= rhs; }
bool BedRange::operator>=(int rhs) const { return range.second >= rhs; }
int BedRange::operator-(int rhs) const { return range.second - rhs; }

Bed::Bed(const std::string &ifile) {
  using namespace RJBUtil;
  Splitter<std::string> files(ifile, ",");
  size_t reserve = 1000000;
  variants_.reserve(reserve);
  for (const auto &f : files) {
    if (!check_file_exists(f)) {
      std::cerr << "No mask file provided or incorrect path to mask file. "
                << f << std::endl;
      return;
    }
    std::ifstream ifs(f);
    std::string line;
    int lineno = -1;

    while (std::getline(ifs, line)) {
      lineno++;
      if (lineno == reserve) {
        reserve += 1000000;
        variants_.reserve(reserve);
      }
      RJBUtil::Splitter<std::string> splitter(line, "\t");
      FileValidator::validate_bed_line(splitter, lineno);

      std::stringstream ss;
      ss << splitter[0] << "," << splitter[1] << "," << splitter[2] << "," << splitter[3] << "," << splitter[4];

      variants_.emplace(ss.str());
    }
  }
}

Bed::Bed(std::stringstream &iss) {
  using namespace RJBUtil;
  std::string line;
  int lineno = -1;

  while (std::getline(iss, line)) {
    lineno++;
    RJBUtil::Splitter<std::string> splitter(line, "\t");
    FileValidator::validate_bed_line(splitter, lineno);

    std::stringstream ss;
    ss << splitter[0] << "," << splitter[1] << "," << splitter[2] << "," << splitter[3] << "," << splitter[4];

    variants_.emplace(ss.str());
  }
}

bool Bed::check_variant(const std::string &variant) {
  return variants_.contains(variant);
}

bool Bed::empty() { return variants_.empty(); }

unsigned long Bed::size() { return variants_.size(); }
