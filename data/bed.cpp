//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include <iostream>

#include "bed.hpp"
#include "../utility/filesystem.hpp"

BedRange::BedRange(std::pair<int, int> &range) : range(range) {}

BedRange::BedRange(int first, int second) : range({first, second}) {}

bool BedRange::operator<(int rhs) const { return range.second < rhs; }
bool BedRange::operator>(int rhs) const { return range.first > rhs; }
bool BedRange::operator==(int rhs) const { return range.first <= rhs && range.second >= rhs; }
bool BedRange::operator<=(int rhs) const { return range.first <= rhs; }
bool BedRange::operator>=(int rhs) const { return range.second >= rhs; }
int BedRange::operator-(int rhs) const { return range.second - rhs; }

Bed::Bed(const std::string &ifile) {
  if(!check_file_exists(ifile)) {
    std::cerr << "No mask file provided or incorrect path."  << std::endl;
    return;
  }
  std::ifstream ifs(ifile);
  std::string line;

  while (std::getline(ifs, line)) {
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	int first = std::stoi(splitter[1]);
	int second = std::stoi(splitter[2]);

	// Initialize vector for chromosome
	if (ranges_.find(splitter[0]) == ranges_.end())
	  ranges_[splitter[0]] = std::vector<BedRange>();

	auto it = std::lower_bound(ranges_[splitter[0]].begin(), ranges_[splitter[0]].end(), first);

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

Bed::Bed(std::stringstream &ss) {
  std::string line;

  while (std::getline(ss, line, '\n')) {
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	int first = std::stoi(splitter[1]);
	int second = std::stoi(splitter[2]);

	// Initialize vector for chromosome
	if (ranges_.find(splitter[0]) == ranges_.end())
	  ranges_[splitter[0]] = std::vector<BedRange>();

	auto it = std::lower_bound(ranges_[splitter[0]].begin(), ranges_[splitter[0]].end(), first);

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

bool Bed::check_variant(const std::string &chr, std::pair<std::string, std::string> &&pos) {
  int start = std::stoi(pos.first);
  int end = std::stoi(pos.second);

  if (ranges_.find(chr) == ranges_.end()) {
    return false;
  }

  auto it = std::lower_bound(ranges_[chr].begin(), ranges_[chr].end(), start); // First range >= start
  if (it == ranges_[chr].end()) {
    return false;
  }
  return *it <= end;
}

bool Bed::check_variant(const std::string &chr, const std::string &pos) {
  int ipos = std::stoi(pos);

  if (ranges_.find(chr) == ranges_.end()) {
    return false;
  }

  auto it = std::lower_bound(ranges_[chr].begin(), ranges_[chr].end(), ipos);
  if (it == ranges_[chr].end()) {
	return false;
  }
  return *it == ipos;
}

bool Bed::check_variant(const std::string &chr, int pos) {
  if (ranges_.find(chr) == ranges_.end()) {
	throw std::out_of_range("Chromosome not in bed file.");
  }

  auto it = std::lower_bound(ranges_[chr].begin(), ranges_[chr].end(), pos);
  if (it == ranges_[chr].end()) {
	return false;
  }
  return *it == pos;
}

bool Bed::empty() {
  return ranges_.empty();
}

unsigned long Bed::size() {
  return ranges_.size();
}

unsigned long Bed::chromosome_count(const std::string &k) {
  return ranges_[k].size();
}
