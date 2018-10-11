//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include "weight.hpp"
#include "../utility/filesystem.hpp"

// TODO: Add check for lines in the wrong format.
Weight::Weight(const std::string &ifile) {
  if (!check_file_exists(ifile)) {
    std::cerr << "No weights provided." << std::endl;
    return;
  }
  std::ifstream ifs(ifile);
  std::string line;

  while (std::getline(ifs, line)) {
    if(line.empty() || line[0] == '#') {
      continue;
    }
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	std::stringstream ss;

	ss << splitter[0] << "-" << splitter[1] << "-" << splitter[2] << "-" << splitter[3];

	double score = std::stod(splitter[4]);

	// Prevent math errors
	scores_[ss.str()] = score;
  }
}

double Weight::get(const std::string &k) const {
  try{
	return scores_.at(k);
  } catch(std::exception &e) {
    std::cerr << "Failed to find weight: " << k << std::endl;
    return 1;
  }
}

double Weight::get(const std::string &k) {
  return scores_.at(k);
}

bool Weight::empty() const {
  return scores_.empty();
}

bool Weight::empty() {
  return scores_.empty();
}

