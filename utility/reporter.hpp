//
// Created by Bohlender,Ryan James on 10/11/18.
//

#ifndef PERMUTE_ASSOCIATE_REPORTER_HPP
#define PERMUTE_ASSOCIATE_REPORTER_HPP

#include <set>

#include "../data/taskqueue.hpp"

class Reporter {
public:
  Reporter(TaskQueue &tq, TaskParams &tp);

  auto report_detail(TaskQueue &tq, TaskParams &tp) -> void;
  auto report_simple(TaskParams &tp) -> void;

private:
  const std::string method_;
  const bool gene_list_;
  const bool testable_;
  std::vector<Result> results_;

  static const std::set<std::string> pvalue_methods_;

  auto extract_results(std::vector<TaskArgs> &tq_results) -> void;
  auto write_to_stream(std::ostream &os, Result &res) -> void;
  auto recalculate_mgit(std::map<std::string, std::map<std::string, Result>> &results) -> void;
};

#endif //PERMUTE_ASSOCIATE_REPORTER_HPP
