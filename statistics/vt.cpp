//
// Created by Bohlender,Ryan James on 10/31/18.
//

#include "vt.hpp"

#include <ctime>

VT_Res::VT_Res() {

}

auto VT_Res::initialize(Gene &gene, arma::vec &pheno, const std::string &ts) -> void{
  arma::sp_mat X(gene.genotypes[ts]);

  arma::uword N = X.n_rows;
  geno_[ts] = arma::vec(X.n_rows * X.n_cols, arma::fill::zeros);
  for(arma::uword i = 0; i < X.n_cols; i++) {
	arma::span span(i * N, std::min((i + 1) * N - 1, geno_[ts].n_elem - 1));
	geno_[ts](span) = arma::vec(X.col(i));
  }

  snpid_[ts] = arma::uvec(X.n_rows * X.n_cols);
  for(arma::uword i = 0; i < X.n_cols; i++) {
	arma::span snp(i * N, (i + 1) * N - 1);
	snpid_[ts](snp).fill(i);
  }

  meanpheno_ = arma::mean(pheno);

  // Values that don't change in permutation
  count_[ts] = arma::vec(X.n_elem, arma::fill::zeros);
  for(arma::uword i = 0; i < X.n_cols; i++) {
	// count_[ts](i) = nCount_[ts](snpid_[ts](i));
    arma::span span(i * N, std::min((i + 1) * N - 1, X.n_elem - 1));
    double val = arma::accu(X.col(i));
    count_[ts](span).fill(val);
  }
  mCount_[ts] = geno_[ts];
  countg_[ts] = arma::vec(count_[ts].n_elem, arma::fill::zeros);
  arma::uword elem = 0;
  arma::vec ucount = arma::unique(count_[ts]);
  for(const auto &v : count_[ts]) {
	countg_[ts](elem) = arma::accu(ucount < v);
	elem++;
  }
  countSquare_[ts] = arma::pow(mCount_[ts], 2);
  arma::vec f = (1 + count_[ts]) / std::sqrt(2 + 2 * N);
  weight_[ts] = 1 / arma::sqrt(f % (1. - f));
  countWeight_[ts] = count_[ts] % weight_[ts];

  // For group
  double nx = countg_[ts].n_elem;
  group_[ts] = arma::vec(nx, arma::fill::ones) + (countg_[ts] - 1);
  ugroup_[ts] = arma::unique(group_[ts]);
  double len = ugroup_[ts].n_elem;
  arr_[ts] = arma::vec(len, arma::fill::zeros);
  oneToLen_[ts] = arma::conv_to<arma::uvec>::from(arma::regspace(0, len-1));
  whiches_[ts] = {};
  for(arma::uword i = 0; i < len; i++) {
	whiches_[ts].push_back(arma::find(group_[ts] == i));
  }

  count_[ts] = sum_groups(mCount_[ts], oneToLen_[ts], ts);
  countSquare_[ts] = sum_groups(countSquare_[ts], oneToLen_[ts], ts);

  // Insensitive to permutation -- TODO move to VT class
  csCount_[ts] = arma::cumsum(count_[ts]);
  csCountSquare_[ts] = arma::cumsum(countSquare_[ts]);
  csCountMeanpheno_[ts] = csCount_[ts] * meanpheno_;
  sqrtCsCountSquare_[ts] = arma::sqrt(csCountSquare_[ts]);

  initialized_[ts] = true;
}

auto VT_Res::is_initialized(const std::string &ts) const -> bool {
  return initialized_.count(ts) > 0;
}

auto VT_Res::get_subset(const std::string &ts) -> arma::uvec & {
  return subset_[ts];
}
auto VT_Res::get_geno(const std::string &ts) -> arma::vec & {
  return geno_[ts];
}

auto VT_Res::sum_groups (arma::vec &X, arma::uvec &range, const std::string &ts) -> arma::vec {
  arma::vec ret(range.n_elem, arma::fill::zeros);
  for(const auto &i : range) {
	ret(i) = arma::sum(X(whiches_[ts][i]));
  }
  return ret;
}

auto VT_Res::get_mCount(const std::string &ts) -> arma::vec & {
  return mCount_[ts];
}

auto VT_Res::get_oneToLen(const std::string &ts) -> arma::uvec & {
  return oneToLen_[ts];
}

auto VT_Res::get_csCountMeanpheno(const std::string &ts) -> arma::vec & {
  return csCountMeanpheno_[ts];
}

auto VT_Res::get_sqrtCsCountSquare(const std::string &ts) -> arma::vec & {
  return sqrtCsCountSquare_[ts];
};

