//
// Created by Bohlender,Ryan James on 10/8/18.
//

#ifndef PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP
#define PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP

#include <vector>
#include <string>
#include <memory>

#include "../data/covariates.hpp"
#include "../data/taskargs.hpp"
#include "../data/taskqueue.hpp"
#include "../data/permutation.hpp"
#include "../data/bed.hpp"
#include "../data/weight.hpp"

class JobDispatcher {
public:
  JobDispatcher(TaskParams &tp);

private:
  // Member functions
  void parse_gene_list();

  // Member variables
  TaskParams tp_;
  TaskQueue tq_;
  Covariates cov_;
  Bed bed_;
  Weight weight_;
  Permute permute_;
  RJBUtil::Splitter<std::string> gene_list_;

  std::shared_ptr<std::vector<std::vector<int32_t>>> permutation_ptr_;
};

#endif //PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP
