#include <catch2/catch.hpp>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <armadillo>

#include "../utility/taskparams.hpp"
#include "../caper/capertask.hpp" // for Stage enum

struct DummyTask {
  long success_threshold;
  long nperm;
  long offset;
};

struct DummyQueue {
  std::vector<DummyTask> tasks;
  size_t nthreads;
  explicit DummyQueue(size_t n) : nthreads(n) {}
  void wait_for_space(size_t) {}
  void dispatch(const DummyTask &t) { tasks.push_back(t); }
  size_t get_nthreads() const { return nthreads; }
};

struct MiniDispatcher {
  TaskParams tp;
  DummyQueue tq;
  Stage stage = Stage::Stage1;
  std::unordered_map<std::string, arma::uword> nvariants;

  explicit MiniDispatcher(TaskParams params)
      : tp(std::move(params)), tq(tp.nthreads - 1) {}

  void multiple_dispatch() {
    if (std::any_of(nvariants.cbegin(), nvariants.cend(),
                    [](const auto &v) { return v.second > 0; })) {
      long total_perm = 0;
      long total_success = 0;
      long perm_step = tp.nperm / tq.get_nthreads();
      long succ_step = tp.success_threshold / tq.get_nthreads();

      for (int i = 0; i < static_cast<int>(tq.get_nthreads()); i++) {
        DummyTask t;
        if (i == static_cast<int>(tq.get_nthreads()) - 1) {
          t = DummyTask{tp.success_threshold - total_success,
                        tp.nperm - total_perm,
                        i * perm_step};
          tq.dispatch(t);
        } else {
          t = DummyTask{succ_step, perm_step, i * perm_step};
          total_perm += perm_step;
          total_success += succ_step;
          tq.dispatch(t);
        }
      }
    }
  }
};

TEST_CASE("multiple_dispatch assigns non-overlapping permutation ranges", "[multiple_dispatch]") {
  TaskParams tp{};
  tp.nperm = 10;
  tp.success_threshold = 6;
  tp.nthreads = 3; // results in 2 worker threads

  MiniDispatcher dispatcher(tp);
  dispatcher.nvariants = {{"tx1", 1}}; // gene has variants

  dispatcher.multiple_dispatch();

  auto &tasks = dispatcher.tq.tasks;
  REQUIRE(tasks.size() == dispatcher.tq.get_nthreads());

  std::sort(tasks.begin(), tasks.end(),
            [](const DummyTask &a, const DummyTask &b) { return a.offset < b.offset; });

  long current = 0;
  long total_perm = 0;
  long total_success = 0;
  for (const auto &t : tasks) {
    REQUIRE(t.offset == current);
    current += t.nperm;
    total_perm += t.nperm;
    total_success += t.success_threshold;
  }

  REQUIRE(current == tp.nperm);
  REQUIRE(total_perm == static_cast<long>(tp.nperm));
  REQUIRE(total_success == static_cast<long>(tp.success_threshold));
}

