//
// Created by Bohlender,Ryan James on 2019-07-30.
//

#ifndef PERMUTE_ASSOCIATE_REFERENCE_HPP
#define PERMUTE_ASSOCIATE_REFERENCE_HPP

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <fstream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <cmath>
#include <array>

/// \brief Container for row(s) from refFlat
struct Gene {
  std::string gene_name;
  std::vector<std::string> transcript;
  std::string chromosome;
  std::string strand;
  std::map<std::string, std::pair<long, long>> positions;
  std::map<std::string, std::pair<long, long>> cds;
  std::map<std::string, long> exons;
  std::map<std::string, std::vector<std::pair<long, long>>> exon_spans;

  bool operator<(Gene &rhs) {
	long left = HUGE_VAL;
	long right = HUGE_VAL;
	for (const auto &v : positions) {
	  if (v.second.first < left) {
		left = v.second.first;
	  }
	}
	for (const auto &v : rhs.positions) {
	  if (v.second.first < right) {
		right = v.second.first;
	  }
	}
	return left < right;
  }
  bool operator<(long pos) {
	long left = HUGE_VAL;
	for (const auto &v : positions) {
	  if (v.second.first < left) {
		left = v.second.first;
	  }
	}
	return left < pos;
  }
  bool operator>(long pos) {
	long left = 0;
	for (const auto &v : positions) {
	  if (v.second.second > left) {
		left = v.second.second;
	  }
	}
	return pos < left ;
  }
  bool operator==(long pos) {
	long left = HUGE_VAL;
	long right = 0;
	for (const auto &v : positions) {
	  if (v.second.first < left) {
		left = v.second.first;
	  }
	}
	for (const auto &v : positions) {
	  if (v.second.second > right) {
		right = v.second.second;
	  }
	}
	return left <= pos && right >= pos;
  }
};

/// \brief Parse and store RefFlat entries
class Reference {
  const int genei = 0;
  const int txi = 1;
  const int chri = 2;
  const int strandi = 3;
  const int starti = 4;
  const int endi = 5;
  const int cds_starti = 6;
  const int cds_endi = 7;
  const int exonsi = 8;
  const int exon_starti = 9;
  const int exon_endi = 10;

  const std::array<std::string, 24> chromosomes {
	"chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
	"chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
	"chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
	"chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
  };
public:
  explicit Reference(const std::string &path);

  // Get genes by position
  std::vector<std::shared_ptr<Gene>> get_gene(const std::string &chr, long pos);

private:
  void parse(std::istream &ifs);

  // Need data structure to be searchable by position, recovering all matching transcripts
  std::map<std::string, std::vector<std::shared_ptr<Gene>>> data_;
  // Searchable by gene name
  std::map<std::string, std::shared_ptr<Gene>> gene_map_;
};

#endif //PERMUTE_ASSOCIATE_REFERENCE_HPP
