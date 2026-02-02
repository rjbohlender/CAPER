//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#ifndef PERMUTE_ASSOCIATE_TASKQUEUE_HPP
#define PERMUTE_ASSOCIATE_TASKQUEUE_HPP

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <random>
#include <thread>
#include <iostream>

#include "taskparams.hpp"

// Three templated arguments to facilitate different run modes
template <typename Operation_t, typename Task_t, typename Reporter_t>
class TaskQueue {

  // Some ideas from
  // https://embeddedartistry.com/blog/2017/2/1/c11-implementing-a-dispatch-queue-using-stdfunction
public:
  std::queue<Operation_t> continue_;
  // Construtors
  TaskQueue(size_t thread_cnt, std::shared_ptr<Reporter_t> reporter,
            TaskParams tp)
      : tp_(tp), threads_(thread_cnt), ntasks_(0), quit_(false),
        nthreads_(thread_cnt), threads_started_(0), reporter_(std::move(reporter)) {
    for (size_t i = 0; i < threads_.size(); i++) {
      threads_[i] = std::thread(std::bind(&TaskQueue::thread_handler, this));
    }
    
    // Wait for all threads to start
    if (nthreads_ > 0) {
        std::unique_lock<std::mutex> lock(lock_);
        cv_started_.wait(lock, [this] { return threads_started_ == (int)nthreads_; });
    }
  }

  // Destructor
  ~TaskQueue() {
    join();
  }

  // Job control
  void join() {
    std::unique_lock<std::mutex> lock(lock_);
    if (quit_) return;

    if (nthreads_ > 0) {
        cv_finish_.wait(lock, [this] { return q_.empty() && ntasks_ == 0; });
    } else {
        // If no worker threads, we must process the queue ourselves or ntasks_ will never reach 0
        while (!q_.empty()) {
            auto op = q_.front();
            q_.pop();
            lock.unlock();
            process_task(op);
            lock.lock();
            ntasks_--;
        }
    }

    quit_ = true;
    cv_.notify_all();
    
    lock.unlock();

    for (size_t i = 0; i < threads_.size(); i++) {
      if (threads_[i].joinable()) {
        threads_[i].join();
      }
    }
  }

  void wait() {
    std::unique_lock<std::mutex> lock(lock_);
    if (quit_) return;

    if (nthreads_ > 0) {
        cv_finish_.wait(lock, [this] { return q_.empty() && ntasks_ == 0; });
    } else {
        while (!q_.empty()) {
            auto op = q_.front();
            q_.pop();
            lock.unlock();
            process_task(op);
            lock.lock();
            ntasks_--;
        }
    }
  }

  // Block until the number of queued tasks is at or below max_size
  void wait_for_space(size_t max_size) {
    std::unique_lock<std::mutex> lock(lock_);
    cv_space_.wait(lock, [this, max_size] { return q_.size() <= max_size; });
  }

  void dispatch(Task_t &args) {
    std::unique_lock<std::mutex> lock(lock_);

    Operation_t op(args, reporter_, rd_(), tp_.verbose);
    q_.push(op);
    ntasks_++;

    lock.unlock();
    cv_.notify_one();
  }

  void dispatch(Task_t &&args) {
    std::unique_lock<std::mutex> lock(lock_);

    Operation_t op(args, reporter_, rd_(), tp_.verbose);
    q_.emplace(op);
    ntasks_++;

    lock.unlock();
    cv_.notify_one();
  }

  void dispatch(Operation_t &op) {
    std::unique_lock<std::mutex> lock(lock_);

    q_.emplace(op);
    ntasks_++;

    lock.unlock();
    cv_.notify_one();
  }

  void dispatch(Operation_t &&op) {
    std::unique_lock<std::mutex> lock(lock_);

    q_.emplace(op);
    ntasks_++;

    lock.unlock();
    cv_.notify_one();
  }

  void unfinished(Operation_t &op) {
    std::lock_guard<std::mutex> lock(cont_lock_);

    continue_.emplace(op);
  }

  void unfinished(Operation_t &&op) {
    std::lock_guard<std::mutex> lock(cont_lock_);

    continue_.emplace(op);
  }

  // Status
  auto empty() -> bool {
    std::lock_guard<std::mutex> lock(lock_);
    return q_.empty();
  }

  auto size() -> size_t {
    std::lock_guard<std::mutex> lock(lock_);
    return q_.size();
  }

  auto get_results() -> std::vector<Task_t> {
    std::lock_guard<std::mutex> lock(lock_);
    return results_;
  }

  auto get_nthreads() const -> const size_t { return nthreads_; }

  auto redispatch() -> void {
    wait();
    std::lock_guard<std::mutex> cont_lock(cont_lock_);
    
    std::unique_lock<std::mutex> lock(lock_);
    while (!continue_.empty()) {
      Operation_t op = continue_.front();
      continue_.pop();
      q_.emplace(op);
      ntasks_++;
    }
    lock.unlock();
    cv_.notify_all();
  }

private:
  // Program arguments
  TaskParams tp_;

  // PRNG
  std::random_device rd_;

  // Threading
  std::mutex lock_;
  std::mutex cont_lock_;
  std::queue<Operation_t> q_;
  std::condition_variable cv_;
  std::condition_variable cv_space_;
  std::condition_variable cv_finish_;
  std::vector<std::thread> threads_;
  bool quit_;
  const size_t nthreads_;

  // Ensuring all jobs are finished
  std::atomic<int> ntasks_;
  std::atomic<int> threads_started_;
  std::condition_variable cv_started_;

  // Messages
  std::shared_ptr<Reporter_t> reporter_;

  // Result storage
  std::vector<Task_t> results_;

  void process_task(Operation_t &op) {
      if (!op.done_) {
        op.run();
        if (!op.done_) {
           unfinished(op);
        } else {
           dispatch(op);
        }
      } else {
        op.finish();
        if (tp_.gene_list) {
          std::lock_guard<std::mutex> lock(lock_);
          results_.emplace_back(op.caperTask);
        }
      }
  }

  auto thread_handler() -> void {
    {
        std::lock_guard<std::mutex> lock(lock_);
        threads_started_++;
    }
    cv_started_.notify_all();

    std::unique_lock lock(lock_);
    while (true) {
      // Wait for data or quit signal
      cv_.wait(lock, [this] { return !q_.empty() || quit_; });

      if (quit_ && q_.empty()) {
          if (ntasks_ == 0) break;
          cv_.wait(lock, [this] { return !q_.empty() || ntasks_ == 0; });
          if (ntasks_ == 0 && q_.empty()) break;
      }

      // After waiting, we have the lock
      if (!q_.empty()) {
        auto op = q_.front();
        q_.pop();
        cv_space_.notify_all();

        // unlock, now that we're done with queue
        lock.unlock();

        process_task(op);

        lock.lock();
        ntasks_--;
        if (ntasks_ == 0 && q_.empty()) {
           cv_finish_.notify_all();
           cv_.notify_all();
        }
      }
    }
  }
};

#endif // PERMUTE_ASSOCIATE_TASKQUEUE_HPP