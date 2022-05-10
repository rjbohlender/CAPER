//
// Created by Bohlender,Ryan James on 8/27/18.
//

#include <algorithm>
#include <set>

#include "covariates.hpp"

#include "../link/binomial.hpp"
#include "../link/gaussian.hpp"
#include "../statistics/bayesianglm.hpp"
#include "../statistics/glm.hpp"
#include "../utility/filevalidator.hpp"
#include "../utility/indexsort.hpp"
#include "../utility/indices.hpp"

Covariates::Covariates(TaskParams tp)
    : tp_(std::move(tp)), nsamples_(0), ncases_(0), ncontrols_(0),
      crand((int)time(nullptr)), linear_(tp_.linear) {
  parse(tp_.covariates_path, tp_.ped_path);
  sorted_ = false;
}

Covariates::Covariates(std::stringstream &ped_ss, std::stringstream &cov_ss)
    : nsamples_(0), ncases_(0), ncontrols_(0), crand((int)time(nullptr)),
      linear_(false) {
  parse(ped_ss, cov_ss);
  sorted_ = false;
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

arma::colvec &Covariates::get_phenotype_vector() { return phenotypes_; }

void Covariates::set_phenotype_vector(const arma::vec &vec) {
  phenotypes_.set_size(vec.n_rows, vec.n_cols);
  phenotypes_ = vec;
}

void Covariates::set_phenotype_vector(const std::vector<int32_t> &vec) {
  phenotypes_ = arma::conv_to<arma::colvec>::from(vec);
}

arma::uword Covariates::get_nsamples() const { return nsamples_; }

arma::uword Covariates::get_ncases() const { return ncases_; }

arma::mat &Covariates::get_covariate_matrix() { return design_; }

arma::vec &Covariates::get_odds() { return odds_; }

arma::vec &Covariates::get_original_phenotypes() { return original_; }

arma::vec &Covariates::get_fitted() { return p_fitted_; }

arma::vec &Covariates::get_mean() { return mean_; }

void Covariates::refit_permuted() {
  bool success;
  if (linear_) {
    Gaussian link("identity");
    GLM<Gaussian> fit(design_, phenotypes_, link, coef_, tp_);
    success = fit.success_;
    p_fitted_ = fit.mu_;
    p_eta_ = fit.eta_;
    p_coef_ = fit.beta_;
    // From Moser and Coombs (2004) -- Get Logistic Regression params without
    // dichotomizing
    double lambda = arma::datum::pi / std::sqrt(3);
    Binomial alt_link("logit");
    arma::vec temp_mu = alt_link.linkinv(
        ((lambda * coef_ / std::sqrt(fit.dev_)) * design_).t());
    p_odds_ =
        temp_mu /
        (1. - temp_mu); // Individual odds, as if the data were dichotomized
  } else {
    Binomial link("logit");
    GLM<Binomial> fit(design_, phenotypes_, link, coef_, tp_);
    success = fit.success_;
    p_odds_ = fit.mu_ / (1. - fit.mu_);
    p_fitted_ = fit.mu_;
    p_eta_ = fit.eta_;
    p_coef_ = fit.beta_;
  }

  // Use BayesianGLM to overcome perfect separation
  if (!success) {
    if (linear_) {
      Gaussian link("identity");
      BayesianGLM<Gaussian> fit(design_, phenotypes_, link);
      p_fitted_ = fit.mu_;
      p_eta_ = fit.eta_;
      p_coef_ = fit.beta_.t();
    } else {
      Binomial link("logit");
      BayesianGLM<Binomial> fit(design_, phenotypes_, link);
      p_odds_ = fit.mu_ / (1. - fit.mu_);
      p_fitted_ = fit.mu_;
      p_eta_ = fit.eta_;
      p_coef_ = fit.beta_.t();
    }
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
  std::ifstream ifs;
  if (!ifile.empty()) {
    ifs = std::ifstream(ifile);
  }
  std::ifstream pfs(pedfile);
  std::string line;

  std::map<std::string, double> sample_phen_map;
  std::vector<double> phenotypes;
  std::vector<std::vector<double>> covariates;
  FileValidator fv;
  int lineno = -1;

  // Parse the phenotype from the ped file.
  while (std::getline(pfs, line)) {
    lineno++;
    if (line[0] == '#') {
      continue;
    }
    RJBUtil::Splitter<std::string> splitter(line, "\t");

    fv.validate_ped_line(splitter, lineno);

    std::string sample_id = splitter[1];
    ped_samples_.insert(sample_id);
    if (ifile.empty()) {
      phenotypes.push_back(std::stod(splitter[5]) - 1);
      if (phenotypes.back() == 1) {
        ncases_++;
      } else {
        ncontrols_++;
      }
    }
    nsamples_++;
    if (linear_) {
      try {
        sample_phen_map[sample_id] = std::stod(splitter[5]);
      } catch (std::exception &e) {
        std::cerr << "Failed to convert quantitative phenotype in .ped file "
                     "column 7.\n";
        throw(e);
      }
    } else {
      try {
        double phen = std::stoi(splitter[5]);
        sample_phen_map[sample_id] = phen - 1;
      } catch (std::exception &e) {
        std::cerr
            << "Failed to convert binary phenotype in .ped file column 6.\n";
        throw(e);
      }
    }
  }
  lineno = -1;
  // Parse the PCA matrix file
  while (ifs.good() && std::getline(ifs, line)) {
    lineno++;
    RJBUtil::Splitter<std::string> splitter(line, " \t");
    fv.validate_cov_line(splitter, lineno);
    if(this->contains(splitter[0])) {
      covariates.emplace_back(std::vector<double>());
      unsigned long i = 0;
      for (const auto &v : splitter) {
        if (i == 0) {
          // Get phenotype of current sample
          auto phen = sample_phen_map[v];

          cov_samples_.push_back(v);

          if (phen == 1) {
            ncases_++;
          } else {
            ncontrols_++;
          }

          phenotypes.push_back(phen);
        } else {
          covariates.back().push_back(std::stod(v));
        }
        i++;
      }
    }
  }
  assert(nsamples_ == ncases_ + ncontrols_ && "Number of samples differs between ped and covariates after parsing.");
  std::cerr << "nsamples: " << nsamples_ << "\n";
  std::cerr << "ncases: " << ncases_ << "\n";
  std::cerr << "ncontrols: " << ncontrols_ << "\n";

  phenotypes_ = arma::conv_to<arma::colvec>::from(phenotypes);
  original_ = arma::conv_to<arma::colvec>::from(phenotypes);

  // Features are j, samples are i
  if (!covariates.empty()) {
    design_ = arma::mat(covariates.size(), covariates[0].size() + 1,
                        arma::fill::zeros);
    for (arma::uword i = 0; i < covariates.size(); i++) {
      for (arma::uword j = 0; j < covariates[0].size() + 1; j++) {
        if (j == 0) {
          // First feature is just 1s, intercept
          // Do this here
          design_(i, j) = 1;
        } else {
          design_(i, j) = covariates[i][j - 1];
        }
      }
    }
  } else {
    design_ = arma::mat(nsamples_, 1, arma::fill::ones);
  }
#ifndef NDEBUG
  std::cerr << "Covariates_ n_rows = " << design_.n_rows
            << " Covariates_ n_cols = " << design_.n_cols << "\n";
#endif
}

/**
 * @brief For debugging
 * @param ss A stringstream to be parsed.
 */
void Covariates::parse(std::stringstream &ped_ss, std::stringstream &cov_ss) {
  std::string line;
  std::map<std::string, double> sample_phen_map;
  std::vector<double> phenotypes;
  std::vector<std::vector<double>> covariates;

  // Parse the phenotype from the ped file.
  while (std::getline(ped_ss, line, '\n')) {
    if (line[0] == '#') {
      continue;
    }
    RJBUtil::Splitter<std::string> splitter(line, "\t");

    std::string sample_id = splitter[1];
    ped_samples_.insert(sample_id);

    nsamples_++;
    if (linear_) {
      try {
        sample_phen_map[sample_id] = std::stod(splitter[5]);
      } catch (std::exception &e) {
        std::cerr << "Failed to convert quantitative phenotype in .ped file "
                     "column 7.\n";
        throw(e);
      }
    } else {
      try {
        double phen = std::stoi(splitter[5]);
        sample_phen_map[sample_id] = phen - 1;
      } catch (std::exception &e) {
        std::cerr
            << "Failed to convert binary phenotype in .ped file column 6.\n";
        throw(e);
      }
    }
  }

  // Parse the PCA matrix file
  while (std::getline(cov_ss, line, '\n')) {
    RJBUtil::Splitter<std::string> splitter(line, " \t");
    if (ped_samples_.count(splitter[0]) > 0) {
      covariates.emplace_back(std::vector<double>());

      unsigned long i = 0;
      for (const auto &v : splitter) {
        if (i == 0) {
          // Get phenotype of current sample
          auto phen = sample_phen_map[v];

          cov_samples_.push_back(v);

          if (phen == 1) {
            ncases_++;
          } else {
            ncontrols_++;
          }

          phenotypes.push_back(phen);
        } else {
          covariates.back().push_back(std::stod(v));
        }
        i++;
      }
    }
  }
  std::cerr << "nsamples: " << nsamples_ << "\n";
  std::cerr << "ncases: " << ncases_ << "\n";
  std::cerr << "ncontrols: " << ncontrols_ << "\n";

  phenotypes_ = arma::conv_to<arma::vec>::from(phenotypes);
  original_ = arma::conv_to<arma::vec>::from(phenotypes);

  // Features are j, samples are i
  if (!covariates.empty()) {
    design_ = arma::mat(covariates.size(), covariates[0].size() + 1,
                        arma::fill::zeros);
    for (arma::uword i = 0; i < covariates.size(); i++) {
      for (arma::uword j = 0; j < covariates[0].size() + 1; j++) {
        if (j == 0) {
          // First feature is just 1s, intercept
          // Do this here
          design_(i, j) = 1;
        } else {
          design_(i, j) = covariates[i][j - 1];
        }
      }
    }
  } else {
    design_ = arma::mat(nsamples_, 1, arma::fill::ones);
  }
#ifndef NDEBUG
  std::cerr << "Covariates_ n_rows = " << design_.n_rows
            << " Covariates_ n_cols = " << design_.n_cols << "\n";
#endif
}

void Covariates::fit_null() {
  if (linear_) {
    Gaussian link("identity");
    GLM<Gaussian> fit(design_, phenotypes_, link, tp_);
    fitted_ = fit.mu_;
    eta_ = fit.eta_;
    coef_ = fit.beta_;
    // From Moser and Coombs (2004) -- Get Logistic Regression params without
    // dichotomizing
    double lambda = arma::datum::pi / std::sqrt(3);
    Binomial alt_link("logit");
    arma::vec temp_mu = alt_link.linkinv(
        ((lambda * coef_ / std::sqrt(fit.dev_)) * design_).t());
    odds_ = temp_mu /
            (1. - temp_mu); // Individual odds, as if the data were dichotomized

    // Initially, permuted values are equal to the first fit.
    p_odds_ = odds_;
    p_fitted_ = fitted_;
    p_eta_ = eta_;
    p_coef_ = coef_;
  } else {
    Binomial link("logit");
    GLM<Binomial> fit(design_, phenotypes_, link, tp_);
    odds_ = fit.mu_ / (1. - fit.mu_);
    fitted_ = fit.mu_;
    eta_ = fit.eta_;
    coef_ = fit.beta_;

    // Initially, permuted values are equal to the first fit.
    p_odds_ = odds_;
    p_fitted_ = fitted_;
    p_eta_ = eta_;
    p_coef_ = coef_;
  }
}

/**
 * @brief Uses the header of the matrix to sort the covariates
 * @param header Header of the matrix file
 *
 * If cov_samples_ is empty, covariates weren't provided, phenotypes are in ped
 * file order. Otherwise, phenotypes are in covariate file order.
 */
void Covariates::sort_covariates(std::string &header) {
  RJBUtil::Splitter<std::string> splitter(header, "\t");

  std::map<std::string, arma::uword> header_map;

  // Cov order
  if (!cov_samples_.empty()) {
    if (splitter.size() - static_cast<int>(Indices::first) !=
        cov_samples_.size()) {
      throw(std::runtime_error("The number of samples differs between the "
                               "covariate file and the matrix file."));
    }
    arma::uvec cov_indices = arma::uvec(cov_samples_.size(), arma::fill::zeros);
    for (auto it = cov_samples_.begin(); it != cov_samples_.end(); it++) {
      header_map[*it] = std::distance(cov_samples_.begin(), it);
    }
    for (arma::uword i = 0; i < phenotypes_.n_rows; i++) {
      cov_indices(i) =
          header_map[splitter[i + static_cast<int>(Indices::first)]];
    }
    header_map.clear();
    arma::uvec ped_indices = arma::uvec(
        ped_samples_.size(), arma::fill::zeros); // Separate ped indices
    for (auto it = ped_samples_.begin(); it != ped_samples_.end(); it++) {
      header_map[*it] = std::distance(ped_samples_.begin(), it);
    }
    for (arma::uword i = 0; i < phenotypes_.n_rows; i++) {
      ped_indices(i) =
          header_map[splitter[i + static_cast<int>(Indices::first)]];
    }
    std::set<std::string> seen;
    for (auto v : cov_samples_) {
      if (seen.find(v) == seen.end()) {
        seen.insert(v);
      } else {
        std::cerr << "cov duplicate: " << v << std::endl;
      }
    }

    seen.clear();
    for (auto v : splitter) {
      if (seen.find(v) == seen.end()) {
        seen.insert(v);
      } else {
        std::cerr << "header duplicate: " << v << std::endl;
      }
    }

    // Sort the phenotypes and covariates according to the order in the matrix
    // file.
    original_ = phenotypes_(cov_indices); // Phenotypes are parsed in covariate
                                          // order if covariates are provided.
    phenotypes_ = phenotypes_(cov_indices);
    design_ = design_.rows(cov_indices);
  } else {
    if (splitter.size() - static_cast<int>(Indices::first) !=
        ped_samples_.size()) {
      throw(std::runtime_error("The number of samples differs between the ped "
                               "file and the matrix file."));
    }
    arma::uvec indices = arma::uvec(ped_samples_.size(), arma::fill::zeros);
    for (auto it = ped_samples_.begin(); it != ped_samples_.end(); it++) {
      header_map[*it] = std::distance(ped_samples_.begin(), it);
    }
    for (arma::uword i = 0; i < phenotypes_.n_rows; i++) {
      indices(i) = header_map[splitter[i + static_cast<int>(Indices::first)]];
    }

    // Sort the phenotypes and covariates according to the order in the matrix
    // file.
    original_ = phenotypes_(indices);
    phenotypes_ = phenotypes_(indices);
    design_ = design_.rows(indices);
  }

  fit_null();

  sorted_ = true;
}

bool Covariates::is_sorted() const { return sorted_; }

arma::vec &Covariates::get_coef() { return p_coef_; }

arma::vec Covariates::get_residuals() const { return phenotypes_ - p_fitted_; }

bool Covariates::contains(const std::string &sample) {
  return (ped_samples_.count(sample) > 0);
}
