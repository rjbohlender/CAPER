//
// Created by Bohlender,Ryan James on 9/20/18.
//

#ifndef PERMUTE_ASSOCIATE_CALCCASM_HPP
#define PERMUTE_ASSOCIATE_CALCCASM_HPP

#include <string>
#include <map>

class CalcCASM {
  double ncase;
  double ncontrol;
public:
  CalcCASM(double ncase, double ncontrol);

  double get_score(const std::string &chr, int spos, int epos, const std::string &type);

private:
  // Hardcoded path to database.
  std::string phastcons_path_ = "../db/phastcons-hg19-vertebrate-cds.txt";

  std::map<std::string, std::map<int, double>> phastcons_scores_;
  static const std::map<std::string, double> cons_controlled_ratio_;

  void parse_phastcons();

};

#endif //PERMUTE_ASSOCIATE_CALCCASM_HPP
