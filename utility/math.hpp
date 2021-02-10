//
// Created by Bohlender,Ryan James on 4/22/20.
//

#ifndef PERMUTE_ASSOCIATE_MATH_HPP
#define PERMUTE_ASSOCIATE_MATH_HPP

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

namespace mp = boost::multiprecision;

/**
 * @brief Log N choose K function using log gamma function
 * @tparam T A floating point type
 * @param n The total number to choose from
 * @param k The number chosen
 * @return T of N choose K on the log scale
 */
template<typename T>
T lcnk(T n, T k) {
  if (n == k) {
	return 0;
  }
  return mp::lgamma(n + 1) - (mp::lgamma(k + 1) + mp::lgamma(n - k + 1));
}

/**
 * @brief Poisson confidence interval for p-values in permutation
 * @tparam T A floating point type
 * @param successes The number of successes counted in permutation
 * @param permutations The total number of permutations
 * @return Pair<T, T> containing the confidence interval at 2.5% and 97.5%
 */
template<typename T>
std::pair<T, T> poisson_ci(T successes, T permutations, double lower = 0.025, double upper = 0.975) {
  T lo_ci;
  T hi_ci;
  boost::math::chi_squared hi_dist(2 * successes + 2);
  if (successes > 0) {
	boost::math::chi_squared lo_dist(2 * successes);
	lo_ci = boost::math::quantile(lo_dist, lower) / 2.;
  } else {
	lo_ci = 0;
  }
  hi_ci = boost::math::quantile(hi_dist, upper) / 2.;

  return std::make_pair(lo_ci / permutations, hi_ci / permutations);
}

/**
 * @brief Calculate the set difference between two iterables
 * @tparam T An iterable
 * @param x An iterable collection of elements, possibly with duplicates
 * @param y An iterable collection of elements, possibly with duplicates
 * @return A modified copy of x
 */
template<typename T>
T setdiff(T x, T y) {
  typename T::iterator it;
  for (size_t j = 0; j < y.size(); j++) {
	while (*(it = std::find(x.begin(), x.end(), y[j])) == y[j]) {
	  x.erase(it);
	}
  }
  return x;
}

/**
 * @brief Unbiased estimator of p in a geometric trial.
 * @tparam T A numeric type.
 * @param m The number of successes.
 * @param n The number of trials.
 * @return The unbiased estimate of p.
 *
 * @note Needed to correct for p estimation in adaptive permutation. Unbiased
 * estimate due to J.B.S. Haldane (Biometrika , Nov., 1945, Vol. 33, No. 3
 * (Nov., 1945), pp. 222-225)
 */
template<typename T>
double geometric_p(T m, T n) {
  return (m - 1.) / (n - 1.);
}

/**
 * @brief Find the percentile of a score given a distribution.
 * @tparam T Floating point numeric type.
 * @param score The score to find the percentile of.
 * @param dist The distribution to find the percentile of the score in.
 * @param greater Should values greater than the score be counted, or values less than the score?
 * @return The percentile of the score.
 */
template<typename T>
double percentile_of_score(T score, std::vector<T> dist, bool greater=false) {
  double ret;
  T init = 0.;
  ret = std::accumulate(dist.begin(), dist.end(), init, [&](T &a, T &b) { if (greater) { return a + (score >= b); } else { return a + (score <= b); } });
  ret = (ret + 1.) / (dist.size() + 1.);
  return ret;
}
#endif //PERMUTE_ASSOCIATE_MATH_HPP
