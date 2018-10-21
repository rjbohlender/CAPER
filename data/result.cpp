//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include "result.hpp"

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
	  odds(NAN) {
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
	  odds(NAN) {
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
	  odds(res.odds) {}

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

  return *this;
}

std::ostream &operator<<(std::ostream &stream, const Result &rhs) {
  stream << std::setw(20) << std::left << rhs.rank;
  stream << std::setw(20) << rhs.gene;
  stream << std::setw(20) << rhs.transcript;
  stream << std::setw(20) << std::setprecision(10) << rhs.original;
  stream << std::setw(20) << std::setprecision(10) << rhs.empirical_p;
  stream << std::setw(20) << std::setprecision(10) << rhs.empirical_midp;
  stream << std::setw(20) << rhs.successes;
  stream << std::setw(20) << rhs.permutations;
  stream << std::setw(20) << rhs.mgit_p;
  stream << std::setw(20) << rhs.mgit_successes;
  stream << std::setw(20) << rhs.odds << std::endl;
  return stream;
}

Result &Result::combine(const Result &res) {
  if(gene != res.gene) {
    throw(std::logic_error("Wrong gene in result combine."));
  }
  if(transcript != res.transcript) {
    throw(std::logic_error("Wrong transcript in result combine."));
  }

  successes += res.successes;
  mid_successes += res.mid_successes;
  permutations += res.permutations;

  // Update empirical p and empirical midp
  // TODO Include randomized permutations
  empirical_p = (1. + successes) / (1. + permutations);
  empirical_midp = (1. + mid_successes) / (1. + permutations);

  // Extend permuted values
  permuted.reserve(permuted.size() + res.permuted.size());
  permuted.insert(permuted.end(), res.permuted.begin(), res.permuted.end());

  return *this;
}

void Result::set_rank(int rank) {
  this->rank = rank;
}

