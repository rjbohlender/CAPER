//
// Created by Bohlender,Ryan James on 9/17/18.
//

#ifndef PERMUTE_ASSOCIATE_MAIN_SUPPORT_HPP
#define PERMUTE_ASSOCIATE_MAIN_SUPPORT_HPP

#include <map>
#include <vector>

#include "../data/result.hpp"
#include "../carva/carvatask.hpp"
#include "taskqueue.hpp"
#include "../data/covariates.hpp"

void calc_mgit_pvalues(std::map<std::string, Result> &results,
					   std::vector<std::string> &transcripts,
					   const std::string &method);

void write_simple(TaskQueue &tq, arma::uword ntranscripts, arma::uword ngenes, TaskParams &tp);

void write_detail(TaskQueue &tq, TaskParams &tp);

#endif //PERMUTE_ASSOCIATE_MAIN_SUPPORT_HPP
