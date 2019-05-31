//
// Created by Bohlender,Ryan James on 9/20/18.
//

#ifndef PERMUTE_ASSOCIATE_CALCCASM_HPP
#define PERMUTE_ASSOCIATE_CALCCASM_HPP

#include <string>
#include <unordered_map>

#include "../utility/taskparams.hpp"

class CalcCASM {
  double ncase;
  double ncontrol;
public:
  CalcCASM(double ncase, double ncontrol, TaskParams &tp);

  auto get_score(const std::string &chr, int spos, int epos, const std::string &type) -> double;

private:
  // Hardcoded path to database.
  std::string phastcons_path_ = "../db/phastcons-hg19-vertebrate-cds.txt";

  std::unordered_map<std::string, std::unordered_map<int, double>> phastcons_scores_;
  static const std::unordered_map<std::string, double> cons_controlled_ratio_;

  auto parse_phastcons(TaskParams &tp) -> void;

  auto calc_snv_score(const std::string &chr, int spos) -> double;
  auto calc_ind_score(const std::string &chr, int spos, int epos) -> double;
};

#endif //PERMUTE_ASSOCIATE_CALCCASM_HPP
