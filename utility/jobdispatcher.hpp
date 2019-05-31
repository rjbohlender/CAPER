//
// Created by Bohlender,Ryan James on 10/8/18.
//

#ifndef PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP
#define PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP

#include <vector>
#include <string>
#include <memory>

#include "reporter.hpp"
#include "../data/covariates.hpp"
#include "../carva/carvatask.hpp"
#include "../data/taskqueue.hpp"
#include "../data/permutation.hpp"
#include "../data/bed.hpp"
#include "../data/weight.hpp"
#include "../data/taskqueue2.hpp"
#include "../carva/carvaop.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

class JobDispatcher {
public:
  explicit JobDispatcher(TaskParams &tp, std::shared_ptr<Reporter> reporter);

private:
  // Member functions
  // Dispatch
  void all_gene_dispatcher(std::istream &gt_stream);
  void single_dispatch(Gene &gene);
  // Gene list
  auto find_gene(const std::string &gene);
  void gene_list_dispatcher(std::istream &gt_stream);
  void multiple_dispatch(Gene &gene);

  // Input parsing
  void new_gene(std::stringstream &ss,
				std::string &line,
				RJBUtil::Splitter<std::string> &split,
				RJBUtil::Splitter<std::string> &vsplit);
  void reset_gene(std::stringstream &ss);
  void add_line(std::stringstream &ss,
				std::string &line,
				RJBUtil::Splitter<std::string> &split,
				RJBUtil::Splitter<std::string> &vsplit);

  // Member variables
  TaskParams tp_;
  TaskQueue2<CARVA_Op, TaskArgs, Reporter> tq_;
  Bed bed_;
  Weight weight_;
  Permute permute_;
  RJBUtil::Splitter<std::string> gene_list_;

  // Gene parsing
  std::ifstream gt_ifs_;
  boost::iostreams::filtering_streambuf<boost::iostreams::input> gt_streambuf;
  std::string header_;
  std::string gene_;
  std::string transcript_;
  Stage stage_;
  std::map<std::string, arma::uword> nvariants_;

  // Counters
  arma::uword ntranscripts_ = 0;
  arma::uword ngenes_ = 0;

  std::shared_ptr<Covariates> cov_;
  std::shared_ptr<std::vector<std::vector<int32_t>>> permutation_ptr_;
};

#endif //PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP
