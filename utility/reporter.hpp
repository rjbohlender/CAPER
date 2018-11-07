//
// Created by Bohlender,Ryan James on 10/11/18.
//

#ifndef PERMUTE_ASSOCIATE_REPORTER_HPP
#define PERMUTE_ASSOCIATE_REPORTER_HPP

#include <set>
#include <mutex>

#include "../data/taskargs.hpp"

class Reporter {
public:
  explicit Reporter(TaskParams &tp);
  Reporter(std::vector<TaskArgs> &res, TaskParams &tp);

  auto report_detail(std::vector<TaskArgs> &res, TaskParams &tp) -> void;
  auto report_simple(TaskParams &tp) -> void;

  auto sync_write_simple(Result &res) -> void;
  auto sync_write_detail(const std::string &d, bool testable) -> void;

private:
  // For thread-safe concurrent writing
  std::mutex lock_;
  std::ofstream simple_file_;
  std::ofstream detail_file_;

  const std::string method_;
  const bool gene_list_;
  const bool testable_;
  std::vector<Result> results_;
  std::vector<std::string> details_;

  static const std::set<std::string> pvalue_methods_;

  auto extract_results(std::vector<TaskArgs> &tq_results, TaskParams &tp) -> void;
  auto write_to_stream(std::ostream &os, Result &res) -> void;
  auto recalculate_mgit(std::map<std::string, std::map<std::string, Result>> &results) -> void;
};

#endif //PERMUTE_ASSOCIATE_REPORTER_HPP
