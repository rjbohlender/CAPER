//
// Created by Bohlender,Ryan James on 7/31/18.
//

#ifndef PERMUTE_ASSOCIATE_COVARIATES_HPP
#define PERMUTE_ASSOCIATE_COVARIATES_HPP

#include "../utility/split.hpp"
#include "../statistics/logistic_regression.hpp"
#include "permutation.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <armadillo>

class Covariates {
public:
  explicit Covariates(const std::string& ifile)
  : nsamples_(0),
    ncases_(0) {

    parse(ifile);
    calculate_odds();

    indices_ = arma::regspace<arma::uvec>(0, nsamples_ - 1);

#if 0
    mean_ = calculate_fisher_mean(static_cast<int32_t>(ncases_), odds_);
    std::cerr << "fisher_mean: " << mean_.t();
    std::cerr << "fisher_mean sum: " << arma::sum(mean_) << "\n";
#endif
  }

  explicit Covariates(std::stringstream& ss)
      : nsamples_(0),
        ncases_(0) {
    parse(ss);
    calculate_odds();

    indices_ = arma::regspace<arma::uvec>(0, nsamples_ - 1);
  }

  Covariates(const Covariates &cov) = default;

  Covariates(Covariates &&cov) noexcept
  : phenotypes_(std::move(cov.phenotypes_)),
    original_(std::move(cov.original_)),
    covariates_(std::move(cov.covariates_)),
    odds_(std::move(cov.odds_)),
    nsamples_(cov.nsamples_),
    ncases_(cov.ncases_)
  {}

  Covariates &operator=(const Covariates &cov) = default;

  void print() {
    for(unsigned long i = 0; i < phenotypes_.n_rows; i++) {
      std::cout << phenotypes_[i];
      for(unsigned long j = 0; j < covariates_.n_cols; j++) {
        std::cout << "\t" << covariates_(i, j);
      }
      std::cout << "\n";
    }
  }

  arma::colvec &get_phenotype_vector() {
    return phenotypes_;
  }

  void set_phenotype_vector(arma::colvec &vec) {
    phenotypes_ = vec;
  }

  void set_phenotype_vector(std::vector<int32_t> &vec) {
    phenotypes_ = arma::conv_to<arma::colvec>::from(vec);
  }

  arma::mat &get_covariate_matrix() {
    return covariates_;
  }

  unsigned long get_nsamples() {
    return nsamples_;
  }

  unsigned long get_ncases() {
    return ncases_;
  }

  arma::colvec &get_odds() {
    return odds_;
  }

  arma::colvec &get_original_phenotypes() {
    return original_;
  }

  arma::vec &get_probability() {
    return prob_;
  }

  arma::uvec &get_indices() {
    return indices_;
  }

  arma::vec &get_mean() {
    return mean_;
  }

  void shuffle() {
    indices_ = arma::shuffle(indices_);
  }

  void shuffle(arma::uvec &ma_carriers) {
    indices_ = arma::regspace<arma::uvec>(0, nsamples_ - 1);
    indices_.rows(ma_carriers) = arma::shuffle(indices_.rows(ma_carriers));
  }

  void clear() {
    phenotypes_.reset();
    original_.reset();
    covariates_.reset();
    odds_.reset();
  }

private:
  arma::colvec phenotypes_;
  arma::colvec original_;
  arma::mat covariates_;
  arma::colvec odds_;
  arma::vec prob_;
  arma::uvec indices_;
  arma::vec mean_; // Mean of MFNCH

  unsigned long nsamples_;
  unsigned long ncases_;

  void parse(const std::string& ifile) {
    std::ifstream ifs(ifile);
    std::string line;

    std::vector<double> phenotypes;
    std::vector<std::vector<double>> covariates;

    while(std::getline(ifs, line)) {
      RJBUtil::Splitter<std::string> splitter(line, "\t");

      unsigned long i = 0;
      for(const auto &v : splitter) {
        if(i == 0) {
          int phen = std::stoi(v);

          if(phen == 1)
            ncases_++;
          nsamples_++;

          phenotypes.push_back(std::stoi(v));
        } else {
          if(covariates.size() < i) {
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
    for(int i = 0; i < covariates.size() + 1; i++) {
      for(int j = 0; j < covariates[0].size(); j++) {
        if(i == 0) {
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
 * @param ifile
 */
  void parse(std::stringstream& ss) {
    std::string line;

    std::vector<double> phenotypes;
    std::vector<std::vector<double>> covariates;

    while(std::getline(ss, line, '\n')) {
      RJBUtil::Splitter<std::string> splitter(line, "\t");

      unsigned long i = 0;
      for(const auto &v : splitter) {
        if(i == 0) {
          int phen = std::stoi(v);

          if(phen == 1)
            ncases_++;
          nsamples_++;

          phenotypes.push_back(std::stoi(v));
        } else {
          if(covariates.size() < i) {
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
    for(int i = 0; i < covariates.size() + 1; i++) {
      for(int j = 0; j < covariates[0].size(); j++) {
        if(i == 0) {
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

  void calculate_odds() {
    LogisticRegression lr(covariates_, phenotypes_);
    odds_ = lr.get_odds();
    prob_ = lr.get_probability();
  }
};

#endif //PERMUTE_ASSOCIATE_COVARIATES_HPP
