//
// Created by Bohlender,Ryan James on 9/3/19.
//

#ifndef PERMUTE_ASSOCIATE_JOINTHREADS_HPP
#define PERMUTE_ASSOCIATE_JOINTHREADS_HPP

#include <vector>
#include <thread>

/**
 * @brief Ensures thread joining during destruction
 *
 * Source -- C++ Concurrency in Action 2nd ed. -- Anthony Williams
 */
class JoinThreads {
  std::vector<std::thread> &threads;
public:
  explicit JoinThreads(std::vector<std::thread> &threads_)
	  : threads(threads_) {}
	  ~JoinThreads() {
    for(unsigned long i = 0; i < threads.size(); ++i) {
      if(threads[i].joinable()) {
        threads[i].join();
      }
    }
  }
};

#endif //PERMUTE_ASSOCIATE_JOINTHREADS_HPP
