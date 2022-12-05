//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#ifndef PERMUTE_ASSOCIATE_CAPEROP_HPP
#define PERMUTE_ASSOCIATE_CAPEROP_HPP

#include "../utility/reporter.hpp"
#include "../utility/taskqueue.hpp"
#include "capertask.hpp"

class CAPEROp {
public:
  CAPEROp(CAPERTask &ct, std::shared_ptr<Reporter> reporter, double seed, bool verbose);
  CAPEROp(CAPERTask &&ct, std::shared_ptr<Reporter> reporter, double seed, bool verbose);

  auto run() -> void;
  auto finish() -> void;

  bool done_;
  bool verbose_;
  CAPERTask carvaTask;
private:

  auto op() -> void;

  static auto check_perm(const TaskParams &tp,
						 double perm_val,
						 long success_threshold,
						 std::pair<const std::string, Result> &v) -> void;

  std::shared_ptr<Reporter> reporter_;
};

#endif // PERMUTE_ASSOCIATE_CAPEROP_HPP
