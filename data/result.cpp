//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include "result.hpp"

Result::Result()
	: gene(""),
	  transcript(""),
	  successes(0),
	  mid_successes(0),
	  mgit_successes(0),
	  permutations(0),
	  min_success_at(0),
	  original(NAN),
	  empirical_p(NAN),
	  empirical_midp(NAN),
	  mgit_p(NAN),
	  done(false) {
  permuted = std::vector<double>();
  // permuted = std::map<double, double>();
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
	  mgit_successes(0) {
  permuted = std::vector<double>();
  // permuted = std::map<double, double>();
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
	  mgit_successes(res.mgit_successes) {}

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

  return *this;
}

std::ostream &operator<<(std::ostream &stream, const Result &rhs) {
  stream << std::setw(20) << rhs.gene;
  stream << std::setw(20) << rhs.transcript;
  stream << std::setw(20) << std::setprecision(10) << rhs.original;
  stream << std::setw(20) << std::setprecision(10) << rhs.empirical_p;
  stream << std::setw(20) << std::setprecision(10) << rhs.empirical_midp;
  stream << std::setw(20) << rhs.successes;
  stream << std::setw(20) << rhs.permutations;
  stream << std::setw(20) << rhs.mgit_p;
  stream << std::setw(20) << rhs.mgit_successes << std::endl;
  return stream;
}

