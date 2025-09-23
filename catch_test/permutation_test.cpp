#include <catch2/catch.hpp>

#include <armadillo>

#include <cstdint>
#include <vector>

#include "../data/permutation.hpp"

TEST_CASE("Permute::build_bins groups odds according to epsilon") {
  arma::colvec odds{0.98, 0.10, 0.92, 0.15};
  arma::colvec odds_copy = odds; // preserve original ordering for mean calculations

  SECTION("multiple bins are created when epsilon is small") {
    double epsilon = 0.1;
    Permute perm(12345);

    perm.build_bins(odds, epsilon);

    REQUIRE(perm.bins_built);
    REQUIRE(perm.nsamples == odds.n_rows);

    arma::uvec expected_sort_idx = arma::sort_index(odds);
    REQUIRE(perm.sort_idx.n_elem == expected_sort_idx.n_elem);
    REQUIRE(arma::all(perm.sort_idx == expected_sort_idx));

    std::vector<int32_t> expected_m{2, 2};
    REQUIRE(perm.m == expected_m);

    std::vector<double> expected_bin_odds{(0.10 + 0.15) / 2.0, (0.92 + 0.98) / 2.0};
    REQUIRE(perm.bin_odds.size() == expected_bin_odds.size());
    for (std::size_t i = 0; i < expected_bin_odds.size(); ++i) {
      REQUIRE(perm.bin_odds[i] == Approx(expected_bin_odds[i]));
    }
  }

  SECTION("large epsilon collapses all odds into a single bin") {
    double epsilon = 1.0;
    Permute perm(12345);

    perm.build_bins(odds_copy, epsilon);

    REQUIRE(perm.bins_built);
    REQUIRE(perm.nsamples == odds_copy.n_rows);
    REQUIRE(perm.sort_idx.n_elem == odds_copy.n_rows);

    std::vector<int32_t> expected_m{static_cast<int32_t>(odds_copy.n_rows)};
    REQUIRE(perm.m == expected_m);
    REQUIRE(perm.bin_odds.size() == 1);
    double expected_average = arma::accu(odds_copy) / static_cast<double>(odds_copy.n_rows);
    REQUIRE(perm.bin_odds.front() == Approx(expected_average));
  }
}
