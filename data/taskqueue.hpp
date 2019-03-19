//
// Created by Bohlender,Ryan James on 8/1/18.
//

#ifndef PERMUTE_ASSOCIATE_TASKQUEUE_HPP
#define PERMUTE_ASSOCIATE_TASKQUEUE_HPP

#include <mutex>
#include <queue>
#include <thread>
#include <chrono>
#include <random>
#include <condition_variable>

#include "taskargs.hpp"
#include "permutation.hpp"
#include "../utility/reporter.hpp"

class TaskQueue {

public:
  // Construtors
  explicit TaskQueue(size_t thread_cnt, std::shared_ptr<Reporter> reporter);
  TaskQueue(size_t thread_cnt, std::shared_ptr<Reporter> reporter, bool verbose);

  // Destructor
  ~TaskQueue();

  // Job control
  void join();
  void dispatch(TaskArgs &ta);
  void dispatch(TaskArgs &&ta);

  // Status
  bool empty();
  size_t size();

  // Detail output
  void duplicate();

  std::vector<TaskArgs> &get_results();
  std::vector<TaskArgs> &get_results_duplicate();
  size_t get_nthreads();

private:
  // PRNG
  std::random_device rd_;
  std::mt19937 gen_;

  // Messages
  bool verbose_;
  std::shared_ptr<Reporter> reporter_;

  // Threading
  std::mutex lock_;
  std::queue<TaskArgs> q_;
  std::condition_variable cv_;
  std::vector<std::thread> threads_;
  bool quit_;
  size_t nthreads_;

  // Ensuring all jobs are finished
  std::atomic<int> ntasks_;

  // Result storage
  std::vector<TaskArgs> results_;
  std::vector<TaskArgs> results_detail_;

  void thread_handler();

  // Support member function
  void check_perm(const TaskParams &tp,
				  double perm_val,
				  int success_threshold,
				  std::pair<const std::string, Result> &v);

  double call_method(Methods &method,
					   Gene &gene,
					   Covariates &cov,
					   arma::vec &phenotypes,
					   TaskParams &tp,
					   const std::string &k,
					   bool shuffle,
					   bool detail);

  void stage_1(TaskArgs &ta);
  void stage_2(TaskArgs &ta);
};

#endif //PERMUTE_ASSOCIATE_TASKQUEUE_HPP
