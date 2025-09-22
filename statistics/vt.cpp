//
// Created by Bohlender,Ryan James on 10/31/18.
//

#include "vt.hpp"

#include <utility>
#include <vector>

VT_Res::VT_Res() {

}

auto VT_Res::initialize(Gene &gene, const arma::vec &pheno,
                        const std::string &ts) -> void {
  const arma::sp_mat &X = gene.genotypes[ts];

  arma::vec allele_counts;
  arma::vec allele_squares;
  if (X.n_cols > 0) {
    allele_counts.zeros(X.n_cols);
    allele_squares.zeros(X.n_cols);
    for (auto it = X.begin(); it != X.end(); ++it) {
      const arma::uword col = it.col();
      const double value = *it;
      allele_counts(col) += value;
      allele_squares(col) += value * value;
    }
  }

  arma::uvec order = allele_counts.n_elem > 0
                         ? arma::sort_index(allele_counts, "ascend")
                         : arma::uvec();
  arma::vec counts_sorted = allele_counts.n_elem > 0 ? allele_counts(order)
                                                     : arma::vec();
  arma::vec squares_sorted = allele_squares.n_elem > 0 ? allele_squares(order)
                                                       : arma::vec();

  std::vector<arma::uword> offsets_vec;
  offsets_vec.reserve(counts_sorted.n_elem);
  for (arma::uword i = 0; i < counts_sorted.n_elem; ++i) {
    if (i + 1 == counts_sorted.n_elem || counts_sorted(i) != counts_sorted(i + 1)) {
      offsets_vec.push_back(i);
    }
  }
  arma::uvec offsets = arma::conv_to<arma::uvec>::from(offsets_vec);

  arma::vec cumulative_counts = arma::cumsum(counts_sorted);
  arma::vec cumulative_squares = arma::cumsum(squares_sorted);

  arma::vec cs_counts_at_thresholds = offsets.n_elem > 0
                                          ? cumulative_counts(offsets)
                                          : arma::vec();
  arma::vec sqrt_cs_count_square = offsets.n_elem > 0
                                        ? arma::sqrt(cumulative_squares(offsets))
                                        : arma::vec();

  double mean_pheno = arma::mean(pheno);
  arma::vec cs_count_mean_pheno = cs_counts_at_thresholds * mean_pheno;

  sorted_indices_[ts] = std::move(order);
  threshold_offsets_[ts] = std::move(offsets);
  csCountMeanpheno_[ts] = std::move(cs_count_mean_pheno);
  sqrtCsCountSquare_[ts] = std::move(sqrt_cs_count_square);

  initialized_[ts] = true;
}

auto VT_Res::is_initialized(const std::string &ts) const -> bool {
  return initialized_.count(ts) > 0;
}

auto VT_Res::accumulate_thresholds(const std::string &ts,
                                   const arma::vec &values) const -> arma::vec {
  const arma::uvec &order = sorted_indices_.at(ts);
  const arma::uvec &offsets = threshold_offsets_.at(ts);

  if (order.n_elem == 0 || values.n_elem == 0) {
    return arma::vec();
  }

  arma::vec sorted_values = values(order);
  arma::vec cumulative = arma::cumsum(sorted_values);
  return cumulative(offsets);
}

auto VT_Res::get_csCountMeanpheno(const std::string &ts) -> arma::vec & {
  return csCountMeanpheno_[ts];
}

auto VT_Res::get_sqrtCsCountSquare(const std::string &ts) -> arma::vec & {
  return sqrtCsCountSquare_[ts];
};

