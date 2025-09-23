#include <catch2/catch.hpp>

#include <armadillo>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <numeric>
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

TEST_CASE("Permute::unpack respects shuffle flag") {
  Permute perm(12345);
  const int successes = 3;
  const int bin_size = 5;

  SECTION("without shuffle produces leading ones followed by zeros") {
    StochasticLib3 rng(777);
    auto result = perm.unpack(successes, bin_size, false, rng);

    std::vector<int8_t> expected{1, 1, 1, 0, 0};
    REQUIRE(result == expected);
    REQUIRE(std::accumulate(result.begin(), result.end(), 0) == successes);
  }

  SECTION("with shuffle retains counts and is deterministic with fixed seed") {
    StochasticLib3 rng(777);
    auto shuffled = perm.unpack(successes, bin_size, true, rng);

    REQUIRE(shuffled.size() == bin_size);
    REQUIRE(std::accumulate(shuffled.begin(), shuffled.end(), 0) == successes);

    auto expected = shuffled;

    StochasticLib3 rng_repeat(777);
    auto reshuffled = perm.unpack(successes, bin_size, true, rng_repeat);

    REQUIRE(reshuffled == expected);
  }
}

TEST_CASE("Permute::fisher_yates matches manual shuffle with same seed") {
  std::vector<int8_t> input{0, 1, 2, 3};
  auto original = input;

  StochasticLib3 rng(777);
  Permute::fisher_yates(input, rng);

  REQUIRE(std::is_permutation(input.begin(), input.end(), original.begin(),
                              original.end()));

  std::vector<int8_t> expected = original;
  StochasticLib3 rng_repeat(777);
  for (int i = static_cast<int>(expected.size()) - 1; i >= 1; --i) {
    auto j = rng_repeat.IRandom(0, i);
    std::swap(expected[i], expected[j]);
  }

  REQUIRE(input == expected);
}

TEST_CASE("Permute generates permutations without rebuilding bins") {
  arma::colvec odds{0.2, 0.4, 0.6, 0.8};
  const double epsilon = 0.01;
  Permute perm(2024);
  perm.build_bins(odds, epsilon);

  REQUIRE(perm.bins_built);
  const auto initial_m = perm.m;
  const auto initial_bin_odds = perm.bin_odds;
  const arma::uvec initial_sort_idx = perm.sort_idx;

  const int nperm = 5;
  const int nthreads = 3;
  const arma::uword nsamples = odds.n_rows;

  auto validate_permutations = [&](const auto &perms,
                                   arma::uword expected_ncases) {
    REQUIRE(perms.size() == static_cast<std::size_t>(nperm));
    for (const auto &p : perms) {
      REQUIRE(p.size() == nsamples);
      const int sum = std::accumulate(p.begin(), p.end(), 0);
      REQUIRE(sum == static_cast<int>(expected_ncases));
    }
  };

  SECTION("ncases is zero") {
    const arma::uword ncases = 0;
    auto permutations =
        std::make_shared<std::vector<std::vector<int8_t>>>();

    perm.generate_permutations(permutations, odds, ncases, nperm, nthreads,
                               epsilon);
    validate_permutations(*permutations, ncases);
    REQUIRE(perm.m == initial_m);
    REQUIRE(perm.bin_odds == initial_bin_odds);
    REQUIRE(arma::all(perm.sort_idx == initial_sort_idx));

    const auto epsilon_result =
        perm.epsilon_permutation(nperm, odds, ncases, "transcript", epsilon);
    validate_permutations(epsilon_result, ncases);
    REQUIRE(perm.m == initial_m);
    REQUIRE(perm.bin_odds == initial_bin_odds);
    REQUIRE(arma::all(perm.sort_idx == initial_sort_idx));

    for (const auto &p : *permutations) {
      REQUIRE(std::all_of(p.begin(), p.end(), [](int8_t value) {
        return value == 0;
      }));
    }
    REQUIRE(epsilon_result == *permutations);
  }

  SECTION("ncases is intermediate") {
    const arma::uword ncases = 2;
    auto permutations =
        std::make_shared<std::vector<std::vector<int8_t>>>();

    perm.generate_permutations(permutations, odds, ncases, nperm, nthreads,
                               epsilon);
    validate_permutations(*permutations, ncases);
    REQUIRE(perm.m == initial_m);
    REQUIRE(perm.bin_odds == initial_bin_odds);
    REQUIRE(arma::all(perm.sort_idx == initial_sort_idx));

    const auto epsilon_result =
        perm.epsilon_permutation(nperm, odds, ncases, "transcript", epsilon);
    validate_permutations(epsilon_result, ncases);
    REQUIRE(perm.m == initial_m);
    REQUIRE(perm.bin_odds == initial_bin_odds);
    REQUIRE(arma::all(perm.sort_idx == initial_sort_idx));
  }

  SECTION("ncases equals total sample size") {
    const arma::uword ncases = nsamples;
    auto permutations =
        std::make_shared<std::vector<std::vector<int8_t>>>();

    perm.generate_permutations(permutations, odds, ncases, nperm, nthreads,
                               epsilon);
    validate_permutations(*permutations, ncases);
    REQUIRE(perm.m == initial_m);
    REQUIRE(perm.bin_odds == initial_bin_odds);
    REQUIRE(arma::all(perm.sort_idx == initial_sort_idx));

    const auto epsilon_result =
        perm.epsilon_permutation(nperm, odds, ncases, "transcript", epsilon);
    validate_permutations(epsilon_result, ncases);
    REQUIRE(perm.m == initial_m);
    REQUIRE(perm.bin_odds == initial_bin_odds);
    REQUIRE(arma::all(perm.sort_idx == initial_sort_idx));

    for (const auto &p : *permutations) {
      REQUIRE(std::all_of(p.begin(), p.end(), [](int8_t value) {
        return value == 1;
      }));
    }
    REQUIRE(epsilon_result == *permutations);
  }
}
