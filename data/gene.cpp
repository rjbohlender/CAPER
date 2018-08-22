//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include <string>

#include "gene.hpp"

Gene::Gene(std::stringstream &ss,
		   unsigned long nsamples,
		   std::map<std::string, unsigned long> &nvariants,
		   const CASM &casm)
	: nsamples_(nsamples),
	  nvariants_(nvariants) {
  parse(ss);
  if (!casm.empty()) {
	for (const auto &k : transcripts_) {
	  weights_[k].reshape(nvariants_[k], 1);
	  arma::uword i = 0;
	  for (const auto &v : positions_[k]) {
		weights_[k](i) = casm.get(v);
		i++;
	  }
	}
  }
}

void Gene::print() {
  for (const auto &v : genotypes_)
	std::cout << v.second;
}

arma::mat &Gene::get_matrix(const std::string &k) {
  return genotypes_[k];
}

std::string &Gene::get_gene() {
  return gene_;
}

std::vector<std::string> &Gene::get_transcripts() {
  return transcripts_;
}

arma::vec &Gene::get_weights(const std::string &k) {
  return weights_[k];
}

unsigned long Gene::get_nvariants(const std::string &k) {
  return nvariants_[k];
}

std::vector<std::string> &Gene::get_positions(const std::string &k) {
  return positions_[k];
}

void Gene::clear() {
  // Set matrix size to 0x0 to free space.
  for (auto &v : genotypes_) {
	v.second.reset();
  }
}

void Gene::parse(std::stringstream &ss) {
  std::string line;
  arma::uword i = 0;

  while (std::getline(ss, line, '\n')) {
	if (i == 0) {
	  header_ = line;
	  i++;
	  continue;
	}
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	if (gene_.empty()) {
	  gene_ = splitter[0];
	}
	auto found = std::find(std::begin(transcripts_), std::end(transcripts_), splitter[1]);
	if (found == std::end(transcripts_)) {
	  // Transcript not found -- add
	  transcripts_.push_back(splitter[1]);
	  // Start with matrix transposed
	  genotypes_[transcripts_.back()] = arma::mat(nsamples_, nvariants_[transcripts_.back()]);
	  weights_[transcripts_.back()] = arma::vec(nvariants_[transcripts_.back()], arma::fill::zeros);
	  // Reset counter on new transcript
	  i = 1;
	}
	if (positions_.find(transcripts_.back()) == positions_.end()) {
	  positions_[transcripts_.back()] = std::vector<std::string>();
	}
	positions_[transcripts_.back()].push_back(splitter[2]);

	for (arma::uword j = 3; j < splitter.size(); j++) {
	  double val = std::stod(splitter[j]);
	  // Handle missing data
	  if (val > 2 || val < 0)
		val = 0;
	  genotypes_[transcripts_.back()](j - 3, i - 1) = val;
	}
	i++;
  }
  // Switch to counting minor allele
  for (auto &v : genotypes_) {
	// For each variant
	for (i = 0; i < v.second.n_cols; i++) {
	  // Check allele frequency
	  if (arma::mean(v.second.col(i)) / 2 > 0.5) {
		for (arma::uword j = 0; j < v.second.n_rows; j++) {
		  switch ((int) v.second(j, i)) {
		  case 0: {
			v.second(j, i) = 2;
			break;
		  }
		  case 1: {
			break;
		  }
		  case 2: {
			v.second(j, i) = 0;
			break;
		  }
		  default: break;
		  }
		}
	  }
	}
  }
}

