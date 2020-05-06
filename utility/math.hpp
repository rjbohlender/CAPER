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
std::pair<T, T> poisson_ci(T successes, T permutations) {
  T lo_ci;
  T hi_ci;
  boost::math::chi_squared hi_dist(2 * successes + 2);
  if (successes > 0) {
	boost::math::chi_squared lo_dist(2 * successes);
	lo_ci = boost::math::quantile(lo_dist, 0.025) / 2.;
  } else {
	lo_ci = 0;
  }
  hi_ci = boost::math::quantile(hi_dist, 0.975) / 2.;

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
#endif //PERMUTE_ASSOCIATE_MATH_HPP
