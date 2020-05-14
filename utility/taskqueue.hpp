//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#ifndef PERMUTE_ASSOCIATE_TASKQUEUE_HPP
#define PERMUTE_ASSOCIATE_TASKQUEUE_HPP

#include <mutex>
#include <queue>
#include <thread>
#include <chrono>
#include <random>
#include <condition_variable>
#include <functional>
#include <atomic>

#include "taskparams.hpp"

// Three templated arguments to facilitate different run modes
template <typename Operation_t, typename Task_t, typename Reporter_t>
class TaskQueue {

	// Some ideas from https://embeddedartistry.com/blog/2017/2/1/c11-implementing-a-dispatch-queue-using-stdfunction
public:
  // Construtors
  TaskQueue(size_t thread_cnt, std::shared_ptr<Reporter_t> reporter, TaskParams tp)
  : tp_(tp),
    threads_(thread_cnt),
    ntasks_(0),
    quit_(false),
    nthreads_(thread_cnt),
    reporter_(std::move(reporter)) {
	for (size_t i = 0; i < threads_.size(); i++) {
	  threads_[i] = std::thread(
		  std::bind(&TaskQueue::thread_handler, this));
	}
  }

  // Destructor
  ~TaskQueue() {
	quit_ = true;
	cv_.notify_all();

	// Wait for threads to finish
	for (size_t i = 0; i < threads_.size(); i++) {
	  if (threads_[i].joinable()) {
		threads_[i].join();
	  }
	}
  }

  // Job control
  void join() {
	using namespace std::chrono_literals;
	// Complete jobs
	while (!q_.empty() || ntasks_ > 0) {
	  std::this_thread::sleep_for(0.1s);
	}

	quit_ = true;
	cv_.notify_all();

	for (size_t i = 0; i < threads_.size(); i++) {
	  if (threads_[i].joinable()) {
		threads_[i].join();
	  }
	}
  }

  void wait() {
	using namespace std::chrono_literals;
	// Complete jobs
	while (!q_.empty() || ntasks_ > 0) {
	  std::this_thread::sleep_for(0.1s);
	}
  }

  void dispatch(Task_t &args) {
	std::unique_lock<std::mutex> lock(lock_);

	Operation_t op(args, reporter_, rd_(), tp_.verbose);
	q_.push(op);
	ntasks_++;

	lock.unlock();
	cv_.notify_all();
  }

  void dispatch(Task_t &&args) {
	std::unique_lock<std::mutex> lock(lock_);

	Operation_t op(args, reporter_, rd_(), tp_.verbose);
	q_.emplace(op);
	ntasks_++;

	lock.unlock();
	cv_.notify_all();
  }

  void dispatch(Operation_t &op) {
	std::unique_lock<std::mutex> lock(lock_);

	q_.emplace(op);
	ntasks_++;

	lock.unlock();
	cv_.notify_all();
  }

  void dispatch(Operation_t &&op) {
	std::unique_lock<std::mutex> lock(lock_);

	q_.emplace(op);
	ntasks_++;

	lock.unlock();
	cv_.notify_all();
  }

  void unfinished(Operation_t &op) {
	std::unique_lock<std::mutex> lock(lock_);

	continue_.emplace(op);

	lock.unlock();
	cv_.notify_all();
  }

  void unfinished(Operation_t &&op) {
	std::unique_lock<std::mutex> lock(lock_);

	continue_.emplace(op);

	lock.unlock();
	cv_.notify_all();
  }

  // Status
  auto empty() -> bool {
	return q_.empty();
  }

  auto size() -> size_t {
    return q_.size();
  }

  auto get_results() -> std::vector<Task_t>& {
	return results_;
  }

  auto get_nthreads() const -> const size_t {
    return nthreads_;
  }

  auto redispatch() -> void {
    using namespace std::chrono_literals;
	// Complete current jobs
	while (!q_.empty() || ntasks_ > 0) {
	  std::this_thread::sleep_for(0.1s);
	}
	std::cerr << "Total jobs remaining: " << continue_.size() << std::endl;
    while(!continue_.empty()) {
      Operation_t op = continue_.front();
      continue_.pop();
      dispatch(op);
    }
  }

private:
  // Program arguments
  TaskParams tp_;

  // PRNG
  std::random_device rd_;

  // Threading
  std::mutex lock_;
  std::queue<Operation_t> q_;
  std::queue<Operation_t> continue_;
  std::condition_variable cv_;
  std::vector<std::thread> threads_;
  bool quit_;
  const size_t nthreads_;

  // Ensuring all jobs are finished
  std::atomic<int> ntasks_;

  // Messages
  std::shared_ptr<Reporter_t> reporter_;

  // Result storage
  std::vector<Task_t> results_;

  auto thread_handler() -> void {
	std::unique_lock<std::mutex> lock(lock_);
	do {
	  // Wait for data
	  cv_.wait(lock, [this] {
		return q_.size() || quit_;
	  });

	  // After waiting, we have the lock
	  if (q_.size() && !quit_) {
		auto op = q_.front();
		q_.pop();

		// unlock, now that we're done with queue
		lock.unlock();

		// Logic for running with OPs
		if (!op.done_) {
		  op.run();
		  if(op.done_) {
		    dispatch(op);
		  } else {
		    unfinished(op);
		  }
		} else {
		  op.finish();

		  if(tp_.gene_list) {
			lock.lock();
			results_.emplace_back(op.ta_);
			lock.unlock();
		  }
		}
		ntasks_--;
		lock.lock();
	  }
	} while (!quit_ || ntasks_ > 0);
  }
};

#endif //PERMUTE_ASSOCIATE_TASKQUEUE_HPP
