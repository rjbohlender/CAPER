//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include "result.hpp"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

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

Result::Result(const std::string &gene, const std::string &transcript)
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
  std::stringstream ci;
  ci << std::defaultfloat << std::setprecision(3) << rhs.empirical_ci.first << "," << std::setprecision(3)
	 << rhs.empirical_ci.second;
  std::stringstream midci;
  midci << std::defaultfloat << std::setprecision(3) << rhs.empirical_midci.first << "," << std::setprecision(3)
		<< rhs.empirical_midci.second;
  stream << std::setw(25) << std::defaultfloat << std::left << rhs.gene << " ";
  stream << std::setw(20) << rhs.transcript;
  stream << std::setw(20) << std::setprecision(8) << rhs.original;
  stream << std::setw(20) << std::setprecision(8) << rhs.exact_p;
  stream << std::setw(20) << std::setprecision(8) << rhs.empirical_p;
  stream << std::setw(20) << ci.str();
  stream << std::setw(20) << std::setprecision(8) << rhs.empirical_midp;
  stream << std::setw(20) << midci.str();
  stream << std::setw(20) << rhs.mgit_p;
  stream << std::setw(20) << rhs.successes;
  stream << std::setw(20) << rhs.mgit_successes;
  stream << std::setw(20) << rhs.permutations;
  if(rhs.output_stats) {
    for(const auto &v : rhs.permuted) {
      stream << std::setw(20) << std::setprecision(5) << v;
    }
  }
  stream << std::endl;
  return stream;
}

Result &Result::combine(const Result &res) {
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
  // TODO Include randomized permutations
  if (std::isfinite(rand_perms)) {
	empirical_p = (1. + successes) / (1. + rand_perms);
	empirical_midp = (1. + mid_successes) / (1. + rand_perms);
  } else {
	empirical_p = (1. + successes) / (1. + permutations);
	empirical_midp = (1. + mid_successes) / (1. + permutations);
  }

  update_ci();
  calc_exact_p();

  // Extend permuted values
  permuted.reserve(permuted.size() + res.permuted.size());
  permuted.insert(permuted.end(), res.permuted.begin(), res.permuted.end());

  return *this;
}

void Result::set_rank(int rank) {
  this->rank = rank;
}

void Result::set_odds(double odds) {
  this->odds = odds;
}

void Result::update_ci() {
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
}

void Result::calc_exact_p(double nmac, double nmaj) {
  double m = permutations;
  double n1 = nmac, n2 = nmaj;
  auto b = static_cast<double>(successes);
  double mt;
  try {
	mt = boost::math::binomial_coefficient<double>(n1 + n2, n1);
  } catch (std::overflow_error &e) {
	exact_p = (b + 1.) / (m + 1.);
	return;
  }

  if (!std::isfinite(mt)) {
	exact_p = (b + 1.) / (m + 1.);
	return;
  }

  if (n1 == n2) {
	mt /= 2;
  }

  auto f = [&](double pt) -> double {
	boost::math::binomial dist(m, pt);
	return boost::math::cdf(dist, b);
  };

  double tol = 1e-8;
  double error_estimate;
  double L1;
  size_t max_refinements = 15;
  double I = boost::math::quadrature::trapezoidal(f, 0., 0.5 / (mt + 1), tol, max_refinements, &error_estimate, &L1);
  exact_p = (b + 1.) / (m + 1.) - I;
}

void Result::calc_exact_p() {
  double m = permutations;
  double n1 = nmac, n2 = nmaj;
  auto mt = boost::math::binomial_coefficient<double>(n1 + n2, n1);
  auto b = static_cast<double>(successes);

  if (n1 == n2) {
	mt /= 2;
  }

  auto f = [&](double pt) -> double {
	boost::math::binomial dist(m, pt);
	return boost::math::cdf(dist, b);
  };

  double tol = 1e-8;
  double error_estimate;
  double L1;
  size_t max_refinements = 15;
  exact_p = boost::math::quadrature::trapezoidal(f, 0., 0.5 / (mt + 1), tol, max_refinements, &error_estimate, &L1);
}



