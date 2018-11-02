//
// Created by Bohlender,Ryan James on 8/27/18.
//

#include <algorithm>

#include "covariates.hpp"

#include "../link/binomial.hpp"
#include "../link/gaussian.hpp"
#include "../statistics/glm.hpp"
#include "../statistics/bayesianglm.hpp"



Covariates::Covariates(const std::string& ifile, const std::string& pedfile, bool linear)
	: nsamples_(0),
	  ncases_(0),
	  crand((int)time(nullptr)),
	  linear_(linear) {

  parse(ifile, pedfile);
  indices_ = arma::regspace<arma::uvec>(0, nsamples_ - 1);
  sorted_ = false;
}

Covariates::Covariates(std::stringstream &ss)
	: nsamples_(0),
	  ncases_(0),
	  crand((int)time(nullptr)),
	  linear_(false) {
  parse(ss);
  fit_null();

  indices_ = arma::regspace<arma::uvec>(0, nsamples_ - 1);
}

void Covariates::print() {
  for (unsigned long i = 0; i < phenotypes_.n_rows; i++) {
	std::cout << phenotypes_[i];
	for (unsigned long j = 0; j < design_.n_cols; j++) {
	  std::cout << "\t" << design_(i, j);
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

arma::uword Covariates::get_nsamples() const {
  return nsamples_;
}

arma::uword Covariates::get_ncases() const {
  return ncases_;
}

arma::mat &Covariates::get_covariate_matrix() {
  return design_;
}

arma::colvec &Covariates::get_odds() {
  return odds_;
}

arma::colvec &Covariates::get_original_phenotypes() {
  return original_;
}

arma::vec &Covariates::get_fitted() {
  return fitted_;
}

arma::uvec &Covariates::get_indices() {
  return indices_;
}

arma::vec &Covariates::get_mean() {
  return mean_;
}

void Covariates::shuffle() {
  // Fisher-Yates shuffle is 3-7x faster than armadillo's shuffle
  for(arma::sword i = indices_.n_elem - 1; i >= 0; i--) {
    arma::sword j = static_cast<arma::sword>(crand.IRandom(0, i));
    arma::uword tmp = indices_[i];
    indices_[i] = indices_[j];
    indices_[j] = tmp;
  }
}

void Covariates::clear() {
  phenotypes_.reset();
  original_.reset();
  design_.reset();
  odds_.reset();
}

/**
 * @brief Parse the PCA_Internal matrix file.
 * @param ifile
 */
void Covariates::parse(const std::string &ifile, const std::string &pedfile) {
  std::ifstream ifs(ifile);
  std::ifstream pfs(pedfile);
  std::string line;

  std::map<std::string, double> sample_phen_map;
  std::vector<double> phenotypes;
  std::vector<std::vector<double>> covariates;

  // Parse the phenotype from the ped file.
  while (std::getline(pfs, line)) {
    if(line[0] == '#') {
      continue;
    }
    RJBUtil::Splitter<std::string> splitter(line, "\t");

    std::string sample_id = splitter[1];
    if(linear_) {
      try {
		sample_phen_map[sample_id] = std::stod(splitter[6]);
      } catch(std::exception &e) {
        std::cerr << "Failed to convert quantitative phenotype in .ped file column 7.\n";
        throw(e);
      }
    } else {
      try {
		double phen = std::stoi(splitter[5]);
		sample_phen_map[sample_id] = (phen == 2) ? 1 : 0;
      } catch(std::exception &e) {
        std::cerr << "Failed to convert binary phenotype in .ped file column 6.\n";
        throw(e);
      }
    }
  }

  // Parse the PCA matrix file
  while (std::getline(ifs, line)) {
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	unsigned long i = 0;
	for (const auto &v : splitter) {
	  if (i == 0) {
	    // Get phenotype of current sample
	    auto phen = sample_phen_map[v];

		samples_.push_back(v);

		if (phen == 1)
		  ncases_++;
		nsamples_++;

		phenotypes.push_back(phen);
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
  design_ = arma::mat(covariates.size() + 1, covariates[0].size());
  for (arma::uword i = 0; i < covariates.size() + 1; i++) {
	for (arma::uword j = 0; j < covariates[0].size(); j++) {
	  if (i == 0) {
		// First feature is just 1s, intercept
		// Do this here
		design_(i, j) = 1;
	  } else {
		design_(i, j) = covariates[i - 1][j];
	  }
	}
  }
#ifndef NDEBUG
  std::cerr << "Covariates_ n_rows = " << design_.n_rows << " Covariates_ n_cols = " << design_.n_cols << "\n";
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
  design_ = arma::mat(covariates.size() + 1, covariates[0].size());
  for (arma::uword i = 0; i < covariates.size() + 1; i++) {
	for (arma::uword j = 0; j < covariates[0].size(); j++) {
	  if (i == 0) {
		// First feature is just 1s, intercept
		// Do this here
		design_(i, j) = 1;
	  } else {
		design_(i, j) = covariates[i - 1][j];
	  }
	}
  }
#ifndef NDEBUG
  std::cerr << "Covariates_ n_rows = " << design_.n_rows << " Covariates_ n_cols = " << design_.n_cols << "\n";
#endif
}

void Covariates::fit_null() {
  if(linear_) {
    Gaussian link("identity");
    GLM<Gaussian> fit(design_, phenotypes_, link);
    fitted_ = fit.mu_;
    eta_ = fit.eta_;
    coef_ = fit.beta_.t();
    // From Moser and Coombs (2004) -- Get Logistic Regression params without dichotomizing
    double lambda = arma::datum::pi / std::sqrt(3);
    Binomial alt_link("logit");
    arma::vec temp_mu = alt_link.linkinv(((lambda * coef_ / std::sqrt(fit.dev_)) * design_).t());
    odds_ = temp_mu / (1. - temp_mu); // Individual odds, as if the data were dichotomized
  } else {
	Binomial link("logit");
	GLM<Binomial> fit(design_, phenotypes_, link);
	odds_ = fit.mu_ / (1. - fit.mu_);
	fitted_ = fit.mu_;
	eta_ = fit.eta_;
	coef_ = fit.beta_.t();
  }
}

/**
 * @brief Uses the header of the matrix to sort the covariates
 * @param header Header of the matrix file
 */
void Covariates::sort_covariates(std::string &header) {
  RJBUtil::Splitter<std::string> splitter(header, "\t");

  arma::uvec indices(phenotypes_.n_rows, arma::fill::zeros);

  arma::uword i = 0;
  for(const auto &v : splitter) {
    if (i >= 3) {
      arma::uword j = 0;
      for(const auto &s : samples_) {
        if(s == v) {
          indices[i - 3] = j;
          break;
        }
        j++;
      }
    }
    i++;
  }

  // Sort the phenotypes and covariates according to the order in the matrix file.
  phenotypes_ = phenotypes_(indices);
  design_ = design_.cols(indices);

  fit_null();

  sorted_ = true;
}

bool Covariates::is_sorted() {
  return sorted_;
}

arma::rowvec &Covariates::get_coef() {
  return coef_;
}

arma::vec Covariates::get_residuals() const {
  return phenotypes_ - fitted_;
}

