//
// Created by Bohlender,Ryan James on 8/27/18.
//

#include "covariates.hpp"

Covariates::Covariates(const std::string &ifile)
	: nsamples_(0),
	  ncases_(0) {

  parse(ifile);
  calculate_odds();

  indices_ = arma::regspace<arma::uvec>(0, nsamples_ - 1);
}

Covariates::Covariates(std::stringstream &ss)
	: nsamples_(0),
	  ncases_(0) {
  parse(ss);
  calculate_odds();

  indices_ = arma::regspace<arma::uvec>(0, nsamples_ - 1);
}

void Covariates::print() {
  for (unsigned long i = 0; i < phenotypes_.n_rows; i++) {
	std::cout << phenotypes_[i];
	for (unsigned long j = 0; j < covariates_.n_cols; j++) {
	  std::cout << "\t" << covariates_(i, j);
	}
	std::cout << "\n";
  }
}

arma::colvec &Covariates::get_phenotype_vector() {
  return phenotypes_;
}

void Covariates::set_phenotype_vector(arma::colvec &vec) {
  phenotypes_ = vec;
}

void Covariates::set_phenotype_vector(std::vector<int32_t> &vec) {
  phenotypes_ = arma::conv_to<arma::colvec>::from(vec);
}

unsigned long Covariates::get_nsamples() {
  return nsamples_;
}

unsigned long Covariates::get_ncases() {
  return ncases_;
}

arma::mat &Covariates::get_covariate_matrix() {
  return covariates_;
}

arma::colvec &Covariates::get_odds() {
  return odds_;
}

arma::colvec &Covariates::get_original_phenotypes() {
  return original_;
}

arma::vec &Covariates::get_probability() {
  return prob_;
}

arma::uvec &Covariates::get_indices() {
  return indices_;
}

arma::vec &Covariates::get_mean() {
  return mean_;
}

void Covariates::shuffle() {
  indices_ = arma::shuffle(indices_);
}

void Covariates::clear() {
  phenotypes_.reset();
  original_.reset();
  covariates_.reset();
  odds_.reset();
}

void Covariates::parse(const std::string &ifile) {
  std::ifstream ifs(ifile);
  std::string line;

  std::vector<double> phenotypes;
  std::vector<std::vector<double>> covariates;

  while (std::getline(ifs, line)) {
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	unsigned long i = 0;
	for (const auto &v : splitter) {
	  if (i == 0) {
		int phen = std::stoi(v);

		if (phen == 1)
		  ncases_++;
		nsamples_++;

		phenotypes.push_back(std::stoi(v));
	  } else {
		if (covariates.size() < i) {
		  covariates.emplace_back(std::vector<double>());
		  covariates.at(i - 1).push_back(std::stod(v));
		} else {
		  covariates.at(i - 1).push_back(std::stod(v));
		}
	  }
	  i++;
	}
  }
  std::cerr << "nsamples_: " << nsamples_ << "\n";
  std::cerr << "ncases: " << ncases_ << "\n";

  phenotypes_ = arma::conv_to<arma::colvec>::from(phenotypes);
  original_ = arma::conv_to<arma::colvec>::from(phenotypes);

  // Features are i, samples are j
  covariates_ = arma::mat(covariates.size() + 1, covariates[0].size());
  for (arma::uword i = 0; i < covariates.size() + 1; i++) {
	for (arma::uword j = 0; j < covariates[0].size(); j++) {
	  if (i == 0) {
		// First feature is just 1s, intercept
		// Do this here
		covariates_(i, j) = 1;
	  } else {
		covariates_(i, j) = covariates[i - 1][j];
	  }
	}
  }
#ifndef NDEBUG
  std::cerr << "Covariates_ n_rows = " << covariates_.n_rows << " Covariates_ n_cols = " << covariates_.n_cols << "\n";
#endif
}

/**
 * @brief For debugging
 * @param ss A stringstream to be parsed.
 */
void Covariates::parse(std::stringstream &ss) {
  std::string line;

  std::vector<double> phenotypes;
  std::vector<std::vector<double>> covariates;

  while (std::getline(ss, line, '\n')) {
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	unsigned long i = 0;
	for (const auto &v : splitter) {
	  if (i == 0) {
		int phen = std::stoi(v);

		if (phen == 1)
		  ncases_++;
		nsamples_++;

		phenotypes.push_back(std::stoi(v));
	  } else {
		if (covariates.size() < i) {
		  covariates.emplace_back(std::vector<double>());
		  covariates.at(i - 1).push_back(std::stod(v));
		} else {
		  covariates.at(i - 1).push_back(std::stod(v));
		}
	  }
	  i++;
	}
  }

#ifndef NDEBUG
  std::cerr << "nsamples_: " << nsamples_ << "\n";
  std::cerr << "ncases_: " << ncases_ << "\n";
#endif

  phenotypes_ = arma::conv_to<arma::colvec>::from(phenotypes);
  original_ = arma::conv_to<arma::colvec>::from(phenotypes);

  // Features are i, samples are j
  covariates_ = arma::mat(covariates.size() + 1, covariates[0].size());
  for (arma::uword i = 0; i < covariates.size() + 1; i++) {
	for (arma::uword j = 0; j < covariates[0].size(); j++) {
	  if (i == 0) {
		// First feature is just 1s, intercept
		// Do this here
		covariates_(i, j) = 1;
	  } else {
		covariates_(i, j) = covariates[i - 1][j];
	  }
	}
  }
#ifndef NDEBUG
  std::cerr << "Covariates_ n_rows = " << covariates_.n_rows << " Covariates_ n_cols = " << covariates_.n_cols << "\n";
#endif
}

void Covariates::calculate_odds() {
  LogisticRegression lr(covariates_, phenotypes_);
  odds_ = lr.get_odds();
  prob_ = lr.get_probability();
  eta_ = lr.get_eta();
}

