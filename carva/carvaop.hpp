//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#ifndef PERMUTE_ASSOCIATE_CARVAOP_HPP
#define PERMUTE_ASSOCIATE_CARVAOP_HPP

#include "carvatask.hpp"
#include "../utility/reporter.hpp"
#include "../utility/taskqueue2.hpp"

class CARVAOp {
public:
  CARVAOp(CARVATask &ta, std::shared_ptr<Reporter> reporter, double seed, bool verbose);

  CARVAOp(const CARVAOp &op);
  CARVAOp(CARVAOp &&op) noexcept;
  CARVAOp &operator=(const CARVAOp &rhs);

  auto run() -> void;
  auto finish() -> void;
  auto is_done() const -> bool;
  auto get_args() -> CARVATask;

private:
  // PRNG
  std::mt19937 gen_;

  auto stage1() -> void;
  auto stage2() -> void;

  auto check_perm(const TaskParams &tp,
				  double perm_val,
				  int success_threshold,
				  std::pair<const std::string, Result> &v) -> void;

  auto call_method(Methods &method,
					 Gene &gene,
					 Covariates &cov,
					 arma::vec &phenotypes,
					 TaskParams &tp,
					 const std::string &k,
					 bool shuffle,
					 bool detail) -> double;


  CARVATask ta_;
  bool done_;
  bool verbose_;
  std::shared_ptr<Reporter> reporter_;
};

#endif //PERMUTE_ASSOCIATE_CARVAOP_HPP
