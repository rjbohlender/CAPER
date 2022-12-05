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
      RJBUtil::Splitter<std::string> splitter(line, "\t");
      FileValidator::validate_bed_line(splitter, lineno);

      std::stringstream ss;
      ss << splitter[0] << "," << splitter[1] << "," << splitter[2] << "," << splitter[3] << "," << splitter[4];

      variants_.emplace(ss.str());
    }
  }
}

Bed::Bed(std::stringstream &ss) {
  using namespace RJBUtil;
  std::string line;
  int lineno = -1;

  while (std::getline(ss, line, '\n')) {
    lineno++;
    Splitter<std::string> splitter(line, "\t");
    FileValidator::validate_bed_line(splitter, lineno);

    int first = std::stoi(splitter[1]);
    int second = std::stoi(splitter[2]);

    // Initialize vector for chromosome
    if (ranges_.find(splitter[0]) == ranges_.end()) {
      ranges_[splitter[0]] = std::vector<BedRange>();
    }

    auto it = std::lower_bound(ranges_[splitter[0]].begin(),
                               ranges_[splitter[0]].end(), first);

    if (it == ranges_[splitter[0]].end()) {
      ranges_[splitter[0]].insert(it, {first, second});
    } else {
      // Check if this is a duplicate
      BedRange br{first, second};
      if ((*it) == first) {
        if (br.range.second > (*it).range.second) {
          (*it).range.second = br.range.second;
        }
      } else {
        ranges_[splitter[0]].emplace(it, std::move(br));
      }
    }
  }
}

bool Bed::check_variant(const std::string &variant) {
  return variants_.count(variant) > 0;
}

bool Bed::empty() { return ranges_.empty(); }

unsigned long Bed::size() { return ranges_.size(); }

unsigned long Bed::chromosome_count(const std::string &k) {
  return ranges_[k].size();
}
