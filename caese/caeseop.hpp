//
// Created by Bohlender,Ryan James on 2019-06-06.
//

#ifndef PERMUTE_ASSOCIATE_CAESEOP_HPP
#define PERMUTE_ASSOCIATE_CAESEOP_HPP

#include "caesetask.hpp"
#include "../utility/reporter.hpp"
#include "../utility/taskqueue.hpp"

class CAESEOp {
public:
  CAESEOp(CAESETask &ct, std::shared_ptr<CAESEReporter> reporter, double seed, bool verbose);

  auto run() -> void;
  auto finish() -> void;
  auto is_done() const -> bool;
  auto get_task() -> CAESETask;

private:
  // PRNG
  std::mt19937 gen_;

  auto effectsize() -> void;

  auto check_perm(const TaskParams &tp,
				  double perm_val,
				  int success_threshold,
				  std::pair<const std::string, Result> &v) -> void;

  auto call_method(CAESETask &ct,
				   arma::vec &phenotypes,
				   TaskParams &tp,
				   const std::string &k) -> double;

  CAESETask ta_;
  bool done_;
  bool verbose_;
  std::shared_ptr<CAESEReporter> reporter_;

};

#endif //PERMUTE_ASSOCIATE_CAESEOP_HPP
