//
// Created by Bohlender,Ryan James on 10/31/18.
//

#include "vt.hpp"

#include <ctime>

VT_Res::VT_Res() {

}

auto VT_Res::initialize(Gene &gene, arma::vec &pheno, const std::string &k) -> void{
  arma::sp_mat X(gene.get_matrix(k));
#if ARMA_VERSION_MAJOR >= 9 && ARMA_VERSION_MINOR >= 6
  geno_[k] = arma::vec(X.as_col());
#else
  geno_[k] = arma::vec(X.n_elem, arma::fill::zeros);
  for(arma::uword i = 0; i < X.n_cols; i++) {
	arma::span span(i * X.n_rows, std::min((i + 1) * X.n_rows - 1, geno.n_elem - 1));
	geno_[k](span) = arma::vec(X.col(i));
  }
#endif
  snpid_[k] = arma::uvec(X.n_elem);
  for(arma::uword i = 0; i < X.n_cols; i++) {
	arma::span snp(i * X.n_rows, (i + 1) * X.n_rows - 1);
	snpid_[k](snp).fill(i);
  }

  meanpheno_ = arma::mean(pheno);
  double N = X.n_rows;

  subset_[k] = arma::find(geno_[k] > 0);
  geno_[k] = geno_[k](subset_[k]);
  snpid_[k] = snpid_[k](subset_[k]);

  // Values that don't change in permutation
  count_[k] = arma::vec(geno_[k].n_elem, arma::fill::zeros);
  nCount_[k] = arma::vec(X.n_cols, arma::fill::zeros);
  for(arma::uword i = 0; i < geno_[k].n_elem; i++) {
	nCount_[k](snpid_[k](i))++;
  }
  for(arma::uword i = 0; i < geno_[k].n_elem; i++) {
	count_[k](i) = nCount_[k](snpid_[k](i));
  }
  mCount_[k] = count_[k];
  countg_[k] = arma::vec(count_[k].n_elem, arma::fill::zeros);
  arma::uword elem = 0;
  for(const auto &v : count_[k]) {
	countg_[k](elem) = arma::accu(arma::unique(count_[k]) < v);
	elem++;
  }
  countSquare_[k] = arma::pow(count_[k], 2);
  arma::vec f = (1 + count_[k]) / std::sqrt(2 + 2 * N);
  weight_[k] = 1 / arma::sqrt(f % (1. - f));
  countWeight_[k] = count_[k] % weight_[k];

  // For group
  double nx = countg_[k].n_elem;
  group_[k] = arma::vec(nx, arma::fill::ones) + (countg_[k] - 1.0);
  ugroup_[k] = arma::unique(group_[k]);
  double len = ugroup_[k].n_elem;
  arr_[k] = arma::vec(len, arma::fill::zeros);
  oneToLen_[k] = arma::conv_to<arma::uvec>::from(arma::regspace(0, len-1));
  whiches_[k] = {};
  for(arma::uword i = 0; i < len; i++) {
	whiches_[k].push_back(arma::find(group_[k] == i));
  }

  count_[k] = sum_groups(count_[k], oneToLen_[k], k);
  countSquare_[k] = sum_groups(countSquare_[k], oneToLen_[k], k);

  // Insensitive to permutation -- TODO move to VT class
  csCount_[k] = arma::cumsum(count_[k]);
  csCountSquare_[k] = arma::cumsum(countSquare_[k]);
  csCountMeanpheno_[k] = csCount_[k] * meanpheno_;
  sqrtCsCountSquare_[k] = arma::sqrt(csCountSquare_[k]);

  initialized_[k] = true;
}

auto VT_Res::is_initialized(const std::string &k) const -> bool {
  return initialized_.count(k) > 0;
}

auto VT_Res::get_subset(const std::string &k) -> arma::uvec & {
  return subset_[k];
}
auto VT_Res::get_geno(const std::string &k) -> arma::vec & {
  return geno_[k];
}

auto VT_Res::sum_groups (arma::vec &X, arma::uvec &range, const std::string &k) -> arma::vec {
  arma::vec ret(range.n_elem, arma::fill::zeros);
  for(const auto &i : range) {
	ret(i) = arma::sum(X(whiches_[k][i]));
  }
  return ret;
}

auto VT_Res::get_mCount(const std::string &k) -> arma::vec & {
  return mCount_[k];
}

auto VT_Res::get_oneToLen(const std::string &k) -> arma::uvec & {
  return oneToLen_[k];
}

auto VT_Res::get_csCountMeanpheno(const std::string &k) -> arma::vec & {
  return csCountMeanpheno_[k];
}

auto VT_Res::get_sqrtCsCountSquare(const std::string &k) -> arma::vec & {
  return sqrtCsCountSquare_[k];
};

