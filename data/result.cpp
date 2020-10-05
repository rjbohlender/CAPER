//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include "result.hpp"
#include "../utility/math.hpp"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/multiprecision/cpp_int.hpp>

namespace mp = boost::multiprecision;

Result::Result()
	: gene(),
	  transcript(),
	  successes(0),
	  mid_successes(0),
	  mgit_successes(0),
	  permutations(0),
	  min_success_at(0),
	  original(NAN),
	  empirical_p(NAN),
	  empirical_midp(NAN),
	  mgit_p(NAN),
	  done(false),
	  skippable(false),
	  permuted(),
	  testable(true),
	  odds(NAN),
	  exact_p(NAN),
	  empirical_ci(NAN, NAN),
	  empirical_midci(NAN, NAN),
	  nmac(NAN),
	  nmaj(NAN),
	  rand_perms(NAN),
	  output_stats(false) {
}

Result::Result(const std::string &gene, const std::string &transcript, bool skippable)
	: gene(gene),
	  transcript(transcript),
	  successes(0),
	  mid_successes(0),
	  permutations(0),
	  min_success_at(-1),
	  original(NAN),
	  empirical_p(NAN),
	  empirical_midp(NAN),
	  mgit_p(NAN),
	  done(false),
	  skippable(skippable),
	  permuted(),
	  mgit_successes(0),
	  testable(true),
	  odds(NAN),
	  exact_p(NAN),
	  empirical_ci(NAN, NAN),
	  empirical_midci(NAN, NAN),
	  nmac(NAN),
	  nmaj(NAN),
	  rand_perms(NAN),
	  output_stats(false) {
}

Result::Result(Result &&res) noexcept
	: gene(std::move(res.gene)),
	  transcript(std::move(res.transcript)),
	  successes(res.successes),
	  mid_successes(res.mid_successes),
	  permutations(res.permutations),
	  min_success_at(res.min_success_at),
	  original(res.original),
	  empirical_p(res.empirical_p),
	  empirical_midp(res.empirical_midp),
	  mgit_p(res.mgit_p),
	  done(res.done),
	  skippable(res.skippable),
	  permuted(std::move(res.permuted)),
	  mgit_successes(res.mgit_successes),
	  testable(res.testable),
	  odds(res.odds),
	  exact_p(res.exact_p),
	  empirical_ci(res.empirical_ci),
	  empirical_midci(res.empirical_midci),
	  nmac(res.nmac),
	  nmaj(res.nmaj),
	  rand_perms(res.rand_perms),
	  output_stats(res.output_stats) {}

Result &Result::operator=(Result &&rhs) noexcept {
  gene = std::move(rhs.gene);
  transcript = std::move(rhs.transcript);
  successes = rhs.successes;
  mid_successes = rhs.mid_successes;
  permutations = rhs.permutations;
  min_success_at = rhs.min_success_at;
  original = rhs.original;
  empirical_p = rhs.empirical_p;
  empirical_midp = rhs.empirical_midp;
  mgit_p = rhs.mgit_p;
  done = rhs.done;
  skippable = rhs.skippable;
  permuted = std::move(rhs.permuted);
  mgit_successes = rhs.mgit_successes;
  testable = rhs.testable;
  odds = rhs.odds;
  exact_p = rhs.exact_p;
  nmac = rhs.nmac;
  nmaj = rhs.nmaj;
  rand_perms = rhs.rand_perms;
  output_stats = rhs.output_stats;

  update_ci();

  return *this;
}

std::ostream &operator<<(std::ostream &stream, Result &rhs) {
  std::cerr << "In result: " << rhs.gene << " " << rhs.transcript << " " << rhs.empirical_p << " " << rhs.successes
			<< " " << rhs.permutations << std::endl;

  std::stringstream ci;
  ci << std::defaultfloat << std::setprecision(3) << rhs.empirical_ci.first << "," << std::setprecision(3)
	 << rhs.empirical_ci.second;
  std::stringstream midci;
  midci << std::defaultfloat << std::setprecision(3) << rhs.empirical_midci.first << "," << std::setprecision(3)
		<< rhs.empirical_midci.second;
  stream << std::setw(25) << std::defaultfloat << std::left << rhs.gene << " ";
  stream << std::setw(20) << rhs.transcript;
  stream << std::setw(30) << std::setprecision(15) << rhs.original;
  // stream << std::setw(20) << std::setprecision(8) << rhs.exact_p;
  stream << std::setw(20) << std::setprecision(8) << rhs.empirical_p;
  stream << std::setw(20) << ci.str();
  stream << std::setw(20) << std::setprecision(8) << rhs.empirical_midp;
  stream << std::setw(20) << midci.str();
  stream << std::setw(20) << rhs.mgit_p;
  stream << std::setw(20) << rhs.successes;
  stream << std::setw(20) << rhs.mgit_successes;
  stream << std::setw(20) << rhs.permutations;
  if (rhs.output_stats) {
	for (const auto &v : rhs.permuted) {
	  stream << std::setw(30) << std::setprecision(15) << v;
	}
  }
  stream << std::endl;
  return stream;
}

Result &Result::combine(const Result &res, const TaskParams &tp) {
  if (gene != res.gene) {
	throw (std::logic_error("Wrong gene in result combine."));
  }
  if (transcript != res.transcript) {
	throw (std::logic_error("Wrong transcript in result combine."));
  }

  successes += res.successes;
  mid_successes += res.mid_successes;
  permutations += res.permutations;
  rand_perms += res.rand_perms;

  // Update empirical p and empirical midp
  if (tp.max_perms) {
	if (permutations < *tp.max_perms) {
	  empirical_p = geometric_p(successes, permutations);
	  empirical_midp = geometric_p(mid_successes, static_cast<double>(permutations));
	} else {
	  empirical_p = (1. + successes) / (1. + permutations);
	  empirical_midp = (1. + mid_successes) / (1. + permutations);
	}
  } else {
	if (permutations < tp.nperm) {
	  empirical_p = geometric_p(successes, permutations);
	  empirical_midp = geometric_p(mid_successes, static_cast<double>(permutations));
	} else {
	  empirical_p = (1. + successes) / (1. + permutations);
	  empirical_midp = (1. + mid_successes) / (1. + permutations);
	}
  }

  update_ci();
  // calc_exact_p();

  // Extend permuted values
  permuted.reserve(permuted.size() + res.permuted.size());
  permuted.insert(permuted.end(), res.permuted.begin(), res.permuted.end());

  return *this;
}

void Result::update_ci() {
#if 0
  // Wilson Score Interval Calculation Code
  double z = 1.96;
  double z2 = z * z;

  double sp = (empirical_p + z2 / (2 * permutations)) / (1. + z2 / permutations);
  double ci = z / (1. + z2 / permutations) * std::sqrt(empirical_p * (1. - empirical_p) / permutations + z2 / std::pow(2 * permutations, 2)); // Wilson score interval
  double ci_low = (sp - ci > 0) ? sp - ci : 0;
  double ci_hi = (sp + ci > 1) ? 1 : sp + ci;
  empirical_ci = std::make_pair(ci_low, ci_hi);

  sp = (empirical_midp + z2 / (2 * permutations)) / (1. + z2 / permutations);
  ci = z / (1. + z2 / permutations) * std::sqrt(empirical_midp * (1. - empirical_midp) / permutations
													+ z2 / std::pow(2 * permutations, 2)); // Wilson score interval
  ci_low = (sp - ci > 0) ? sp - ci : 0;
  ci_hi = (sp + ci > 1) ? 1 : sp + ci;
  empirical_midci = std::make_pair(ci_low, ci_hi);
#endif
  empirical_ci = poisson_ci(static_cast<double>(successes), static_cast<double>(permutations));
  empirical_midci = poisson_ci(static_cast<double>(mid_successes), static_cast<double>(permutations));
}

/**
 * @brief Calculate the exact p-value following Phipson and Smyth
 * @param n1
 * @param n2
 * @note Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
 * 		 Authors: Phipson B, Smyth GK., 2010
 */
void Result::calc_exact_p(double n1, double n2) {
  return;
  mp::cpp_bin_float_100 m = permutations;
  mp::cpp_bin_float_100 b = successes;
  mp::cpp_bin_float_100 mt = boost::math::binomial_coefficient<mp::cpp_bin_float_100>(n1 + n2, n1);
#if 0
  if (!std::isfinite(mt)) {
	exact_p = (b + 1.) / (m + 1.);
	return;
  }
#endif

  if (n1 == n2) {
	mt /= 2;
  }

  auto f = [&](mp::cpp_bin_float_100 pt) -> mp::cpp_bin_float_100 {
	boost::math::binomial_distribution<mp::cpp_bin_float_100> dist(m, pt);
	return boost::math::cdf(dist, b);
  };

  mp::cpp_bin_float_100 tol = 1e-8;
  mp::cpp_bin_float_100 error_estimate;
  mp::cpp_bin_float_100 L1;
  size_t max_refinements = 15;
  mp::cpp_bin_float_100 I = boost::math::quadrature::trapezoidal(f,
																 mp::cpp_bin_float_100(0.),
																 mp::cpp_bin_float_100(0.5)
																	 / mp::cpp_bin_float_100(mt + 1),
																 tol,
																 max_refinements,
																 &error_estimate,
																 &L1);
  exact_p = double(mp::cpp_bin_float_100(b + 1.) / mp::cpp_bin_float_100(m + 1.) - I);
}

/**
 * @brief Calculate the exact p-value following Phipson and Smyth
 * @note Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
 * 		 Authors: Phipson B, Smyth GK., 2010
 */
void Result::calc_exact_p() {
  return;
  mp::cpp_bin_float_100 m = permutations;
  double n1 = nmac, n2 = nmaj;
  mp::cpp_bin_float_100 mt = boost::math::binomial_coefficient<mp::cpp_bin_float_100>(n1 + n2, n1);
  mp::cpp_bin_float_100 b = successes;

  if (n1 == n2) {
	mt /= 2;
  }

  auto f = [&](mp::cpp_bin_float_100 pt) -> mp::cpp_bin_float_100 {
	boost::math::binomial_distribution<mp::cpp_bin_float_100> dist(m, pt);
	return boost::math::cdf(dist, b);
  };

  mp::cpp_bin_float_100 tol = 1e-8;
  mp::cpp_bin_float_100 error_estimate;
  mp::cpp_bin_float_100 L1;
  size_t max_refinements = 15;
  exact_p = double(boost::math::quadrature::trapezoidal(f,
														mp::cpp_bin_float_100(0.),
														mp::cpp_bin_float_100(0.5) / mp::cpp_bin_float_100(mt + 1),
														tol,
														max_refinements,
														&error_estimate,
														&L1));
}



