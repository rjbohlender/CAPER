//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include "casm.hpp"

CASM::CASM(const std::string &ifile) {
  std::ifstream ifs(ifile);
  std::string line;

  while (std::getline(ifs, line)) {
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	std::stringstream ss;

	ss << splitter[0] << "-" << splitter[1] << "-" << splitter[2] << "-" << splitter[3];

	double score = std::stod(splitter[4]);

	// Prevent math errors
	scores_[ss.str()] = score > 0 ? std::log(score) : score;
  }
}

CASM::CASM(std::stringstream &iss, bool log_vals) {
  std::string line;

  while (std::getline(iss, line, '\n')) {
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	std::stringstream ss;

	ss << splitter[0] << "-" << splitter[1] << "-" << splitter[2] << "-" << splitter[3];

	double score = std::stod(splitter[4]);

	if (log_vals) {
	  // Prevent math errors
	  scores_[ss.str()] = score > 0 ? std::log(score) : score;
	} else {
	  scores_[ss.str()] = score;
	}
  }
}

double CASM::get(const std::string &k) const {
  return scores_.at(k);
}

double CASM::get(const std::string &k) {
  return scores_.at(k);
}

bool CASM::empty() const {
  return scores_.empty();
}

bool CASM::empty() {
  return scores_.empty();
}

