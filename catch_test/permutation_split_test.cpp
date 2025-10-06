#include <catch2/catch.hpp>
#define private public
#include "../caper/caperop.hpp"
#undef private
#include "../data/result.hpp"
#include "../utility/taskparams.hpp"

TEST_CASE("check_done respects task termination with gene list") {
  TaskParams tp{};
  tp.nperm = 11;   // prime number of permutations
  tp.nthreads = 5; // odd, results in 4 worker threads
  tp.gene_list = std::string("GeneA");

  Result res("GeneA", "Transcript1", false);

  // With fewer permutations than assigned, should not be done
  res.permutations = 4; // termination for last thread will be 5
  CAPEROp::check_done(tp, 1, res, 5);
  REQUIRE_FALSE(res.done);

  // Once permutations reach termination, result should be done
  res.permutations = 5;
  res.done = false;
  CAPEROp::check_done(tp, 1, res, 5);
  REQUIRE(res.done);
}

TEST_CASE("permutation splitting totals match requested count") {
  TaskParams tp{};
  std::vector<std::pair<unsigned long, unsigned int>> configs{
      {11, 3}, // prime permutations, odd threads
      {17, 5}, // prime permutations, odd threads
      {9, 3}   // odd permutations, odd threads
  };
  for (auto [nperm, nthreads] : configs) {
    tp.nperm = nperm;
    tp.nthreads = nthreads;
    tp.gene_list = std::string("GeneA");

    unsigned int workers = tp.nthreads - 1;
    unsigned long base = tp.nperm / workers;
    unsigned long remainder = tp.nperm - base * workers;
    std::vector<unsigned long> parts(workers, base);
    parts.back() += remainder;

    unsigned long total = 0;
    for (auto part : parts) {
      Result r("GeneA", "Transcript1", false);
      r.permutations = part;
      CAPEROp::check_done(tp, 1, r, part);
      REQUIRE(r.done);
      total += r.permutations;
    }
    REQUIRE(total == tp.nperm);
  }
}

TEST_CASE("check_done handles single-thread gene list with max perms") {
  TaskParams tp{};
  tp.nthreads = 1;
  tp.gene_list = std::string("GeneA");
  tp.max_perms = arma::uword(12);

  Result res("GeneA", "Transcript1", false);

  SECTION("below termination threshold remains unfinished") {
    res.permutations = 5;
    CAPEROp::check_done(tp, 1, res, 6);
    REQUIRE_FALSE(res.done);
  }

  SECTION("reaching termination completes result") {
    res.permutations = 6;
    CAPEROp::check_done(tp, 1, res, 6);
    REQUIRE(res.done);
  }
}
