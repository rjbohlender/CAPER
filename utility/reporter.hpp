//
// Created by Bohlender,Ryan James on 10/11/18.
//

#ifndef PERMUTE_ASSOCIATE_REPORTER_HPP
#define PERMUTE_ASSOCIATE_REPORTER_HPP

#include <set>
#include <mutex>

#include "../caper/capertask.hpp"
#include "../power/powertask.hpp"
#include "../caese/caesetask.hpp"

struct ResultLine {
  std::string gene;
  std::string transcript;
  long double original;
  // double exact_p;
  long double empirical_p;
  std::pair<double, double> empirical_p_ci;
  long double empirical_midp;
  std::pair<double, double> empirical_midp_ci;
  double mgit;
  double mgit_midp;
  unsigned long successes;
  unsigned long mgit_successes;
  unsigned long mgit_midp_successes;
  unsigned long permutations;
  std::vector<std::string> stats;
};

class Reporter {
public:
  explicit Reporter(TaskParams &tp);
  Reporter(std::vector<CAPERTask> &res, TaskParams &tp);

  auto report(std::vector<CAPERTask> &&res, TaskParams &tp) -> void;
  auto report_detail(std::vector<CAPERTask> &res, TaskParams &tp) -> void;
  auto report_simple(TaskParams &tp) -> void;
  auto report_vaast(std::vector<CAPERTask> &res, TaskParams &tp) -> void;
  void vaast_sample_index_map(const std::vector<CAPERTask> &res);
  void vaast_sample_index_map(std::vector<std::string> &&samples);

  auto cleanup(TaskParams &tp) -> void;

  auto sync_write_simple(std::unordered_map<std::string, Result> &results, const TaskParams &tp) -> void;
  auto sync_write_detail(const std::string &d, bool gene_testable) -> void;
  auto sync_write_vaast(CAPERTask &ct, const TaskParams &tp) -> void;

  auto sort_simple(const TaskParams &tp) -> void;

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
  const bool pvalues_;
  const bool gene_list_;
  const bool print_testable_;
  std::vector<std::string> details_;
  std::map<std::string, std::string> vaast_;
  std::map<std::string, std::map<std::string, Result>> results_;
  // Holds unfinished genes
  std::vector<std::string> unfinished_;

  static const std::set<std::string> pvalue_methods_;

  auto extract_results(std::vector<CAPERTask> &tq_results, TaskParams &tp) -> void;
  auto write_to_stream(std::ostream &os, Result &res) -> void;
  auto recalculate_mgit(std::map<std::string, std::map<std::string, Result>> &results) -> void;
  auto recalculate_mgit(std::unordered_map<std::string, Result> &results) -> void;
};

class PowerReporter {
public:
  explicit PowerReporter(TaskParams &tp);
  PowerReporter(std::vector<PowerTask> &res, TaskParams &tp);

  auto report(std::vector<PowerTask> &resv, TaskParams &tp) -> void;
  auto report_power(std::vector<PowerTask> &resv, TaskParams &tp) -> void;

  auto cleanup(TaskParams &tp) -> void;

  auto sync_write_power(std::vector<PowerRes> &prv) -> void;

  auto set_ncases(int ncases) -> void;
  auto set_ncontrols(int ncontrols) -> void;

private:
  // For thread-safe concurrent writing
  std::mutex lock_;

  int ncases_;
  int ncontrols_;

  std::ofstream power_file_;

  const std::string method_;

  std::vector<PowerRes> results_;
  std::vector<std::string> details_;
  // Holds unfinished genes
  std::vector<std::string> unfinished_;
};

class CAESEReporter {
public:
  explicit CAESEReporter(TaskParams &tp);
  CAESEReporter(std::vector<CAESETask> &res, TaskParams &tp);

  auto report(std::vector<CAESETask> &resv, TaskParams &tp) -> void;
  auto report_caese(std::vector<CAESETask> &resv, TaskParams &tp) -> void;

  auto cleanup(TaskParams &tp) -> void;

  auto sync_write_caese(std::map<std::string, Result> &prv) -> void;

  auto set_ncases(int ncases) -> void;
  auto set_ncontrols(int ncontrols) -> void;

private:
  // For thread-safe concurrent writing
  std::mutex lock_;

  int ncases_;
  int ncontrols_;

  std::ofstream caese_file_;

  const std::string method_;

  std::vector<Result> results_;
  std::vector<std::string> details_;
  // Holds unfinished genes
  std::vector<std::string> unfinished_;
};
#endif //PERMUTE_ASSOCIATE_REPORTER_HPP
