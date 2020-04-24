//
// Created by Bohlender,Ryan James on 4/22/20.
//

#ifndef PERMUTE_ASSOCIATE_MATH_HPP
#define PERMUTE_ASSOCIATE_MATH_HPP

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

namespace mp = boost::multiprecision;

template<typename T>
T lcnk(T n, T k) {
  if (n == k) {
    return 0;
  }
  return mp::lgamma(n + 1) - (mp::lgamma(k + 1) + mp::lgamma(n - k + 1));
}

template<typename T>
std::pair<T, T> ci(T n, T k, T successes, T permutations) {
  typedef mp::number<mp::cpp_bin_float<1000, mp::digit_base_10, std::allocator<void> > > cpp_bin_float_1000;
  cpp_bin_float_1000 n_ = n;
  cpp_bin_float_1000 k_ = k;
  cpp_bin_float_1000 log_comb = lcnk(n_, k_);

  cpp_bin_float_1000 lo_ci;
  cpp_bin_float_1000 hi_ci;
  boost::math::chi_squared hi_dist(2 * successes + 2);
  if(successes > 0) {
	boost::math::chi_squared lo_dist(2 * successes);
	lo_ci = boost::math::quantile(lo_dist, 0.025) / 2.;
  } else {
	lo_ci = 0;
  }
  hi_ci = boost::math::quantile(hi_dist, 0.975) / 2.;

  // 100% guaranteed underflows for any reasonable sample size using double
  cpp_bin_float_1000 prob = 1. / exp(log_comb);

  return std::make_pair(static_cast<double>(prob + (1. - prob) * (lo_ci / permutations)), static_cast<double>(prob + (1. - prob) * (hi_ci / permutations)));
}

#endif //PERMUTE_ASSOCIATE_MATH_HPP
