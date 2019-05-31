//
// Created by Bohlender,Ryan James on 10/11/18.
//

#ifndef PERMUTE_ASSOCIATE_REPORTER_HPP
#define PERMUTE_ASSOCIATE_REPORTER_HPP

#include <set>
#include <mutex>

#include "../carva/carvatask.hpp"
#include "../statistics/power.hpp"

struct ResultLine {
  std::string gene;
  std::string transcript;
  double original;
  double empirical_p;
  double empirical_midp;
  double mgit;
  unsigned long successes;
  unsigned long mgit_successes;
  unsigned long permutations;
};

class Reporter {
public:
  explicit Reporter(TaskParams &tp);
  Reporter(std::vector<TaskArgs> &res, TaskParams &tp);

  auto report(std::vector<TaskArgs> &res, TaskParams &tp) -> void;
  auto report_detail(std::vector<TaskArgs> &res, TaskParams &tp) -> void;
  auto report_simple(TaskParams &tp) -> void;
  auto report_power(std::vector<TaskArgs> &resv, TaskParams &tp) -> void;
  auto report_vaast(std::vector<TaskArgs> &res, TaskParams &tp) -> void;

  auto sync_write_simple(std::unordered_map<std::string, Result> &results, TaskParams &tp, bool top_only) -> void;
  auto sync_write_detail(const std::string &d, bool testable) -> void;
  auto sync_write_power(std::vector<PowerRes> &prv) -> void;
  auto sync_write_vaast(std::unordered_map<std::string, Result> &results, TaskParams &tp) -> void;

  auto sort_simple(TaskParams &tp) -> void;

  auto set_ncases(int ncases) -> void;
  auto set_ncontrols(int ncontrols) -> void;

private:
  // For thread-safe concurrent writing
  std::mutex lock_;
  // For final simple output
  int ncases_;
  int ncontrols_;
  // Paths -- .tmp for sorting
  std::stringstream simple_path_ss;
  std::stringstream simple_path_tmp_ss;
  std::stringstream detail_path_ss;
  std::stringstream vaast_path_ss;

  std::ofstream simple_file_tmp_;
  std::ofstream simple_file_;
  std::ofstream detail_file_;
  std::ofstream power_file_;
  std::ofstream vaast_file_;

  const std::string method_;
  const bool gene_list_;
  const bool testable_;
  std::vector<Result> results_;
  std::vector<std::string> details_;
  // Holds unfinished genes
  std::vector<std::string> unfinished_;

  static const std::set<std::string> pvalue_methods_;

  auto extract_results(std::vector<TaskArgs> &tq_results, TaskParams &tp) -> void;
  auto write_to_stream(std::ostream &os, Result &res) -> void;
  auto recalculate_mgit(std::map<std::string, std::map<std::string, Result>> &results) -> void;
  auto recalculate_mgit(std::unordered_map<std::string, Result> &results) -> void;
};

#endif //PERMUTE_ASSOCIATE_REPORTER_HPP
