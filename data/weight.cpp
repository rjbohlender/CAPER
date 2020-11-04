//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include <boost/algorithm/string/predicate.hpp>
#include "weight.hpp"
#include "../utility/filesystem.hpp"
#include "../utility/filevalidator.hpp"

Weight::Weight(const std::string &ifile) {
  if (!check_file_exists(ifile)) {
    std::cerr << "No weights provided." << std::endl;
    return;
  }
  std::ifstream ifs(ifile);
  std::string line;
  int lineno = -1;
  FileValidator fileValidator;

  while (std::getline(ifs, line)) {
	lineno++;
    if(line.empty() || boost::starts_with(line, "#")) {
      continue;
    }

	RJBUtil::Splitter<std::string> splitter(line, "\t");
	fileValidator.validate_weight_line(splitter, lineno);

	if(splitter.size() < 5) {
	  std::cerr << "Incorrectly formatted weight line. Line was: " << line << std::endl;
	  std::cerr << "Line should be tab separated and formatted as <chrom> <start_pos> <end_pos> <type> <weight>" << std::endl;
	  throw(std::runtime_error("Incorrect line in weight file."));
    }

	std::stringstream ss;

    if(splitter[3] == "SNP") {
	  ss << splitter[0] << "-" << splitter[1] << "-" << splitter[2] << "-" << "SNV";
    } else {
	  ss << splitter[0] << "-" << splitter[1] << "-" << splitter[2] << "-" << splitter[3];
    }


	double score;
	try {
		score = std::stod(splitter[4]);
	} catch(std::exception &e) {
		std::cerr << "Failed to convert weight to double. Line was: " << line << std::endl;
		std::cerr << "Line should be tab separated and formatted as <chrom> <start_pos> <end_pos> <type> <weight>" << std::endl;
		throw(e);
	}

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

Weight::Weight(std::stringstream &ifile) {
  std::string line;

  while (std::getline(ifile, line, '\n')) {
	if(line.empty() || line[0] == '#') {
	  continue;
	}
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	if(splitter.size() < 5) {
	  std::cerr << "Incorrectly formatted weight line. Line was: " << line << std::endl;
	  std::cerr << "Line should be tab separated and formatted as <chrom> <start_pos> <end_pos> <type> <weight>" << std::endl;
	  throw(std::runtime_error("Incorrect line in weight file."));
	}

	std::stringstream ss;

	if(splitter[3] == "SNP") {
	  ss << splitter[0] << "-" << splitter[1] << "-" << splitter[2] << "-" << "SNV";
	} else {
	  ss << splitter[0] << "-" << splitter[1] << "-" << splitter[2] << "-" << splitter[3];
	}


	double score;
	try {
	  score = std::stod(splitter[4]);
	} catch(std::exception &e) {
	  std::cerr << "Failed to convert weight to double. Line was: " << line << std::endl;
	  std::cerr << "Line should be tab separated and formatted as <chrom> <start_pos> <end_pos> <type> <weight>" << std::endl;
	  throw(e);
	}

	// Prevent math errors
	scores_[ss.str()] = score;
  }
}

