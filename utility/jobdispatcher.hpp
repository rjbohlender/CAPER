//
// Created by Bohlender,Ryan James on 10/8/18.
//

#ifndef PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP
#define PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP

#include <vector>
#include <string>
#include <memory>
#include <chrono>

#include "reporter.hpp"
#include "../data/covariates.hpp"
#include "../carva/carvatask.hpp"
#include "taskqueue.hpp"
#include "../data/permutation.hpp"
#include "../data/bed.hpp"
#include "../data/weight.hpp"
#include "taskqueue2.hpp"
#include "../carva/carvaop.hpp"

#include "main_support.hpp"
#include "reporter.hpp"
#include "filesystem.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

template<typename Operation_t, typename Task_t, typename Reporter_t>
class JobDispatcher {
public:
  JobDispatcher(TaskParams &tp, std::shared_ptr<Reporter_t> reporter)
	  : tp_(tp), tq_(tp_.nthreads - 1, reporter, tp),
		cov_(std::make_shared<Covariates>(tp.covariates_path, tp.ped_path, tp.linear)) {
	// Initialize bed and weights
	if (tp.gene_list)
	  gene_list_ = RJBUtil::Splitter<std::string>(*tp.gene_list, ",");
	if (tp.bed)
	  bed_ = Bed(*tp.bed);
	if (tp.weight)
	  weight_ = Weight(*tp.weight);

	// Update case/control count for reporter
	if(!tp.power) {
	  reporter->set_ncases(cov_->get_ncases());
	  reporter->set_ncontrols(cov_->get_nsamples() - cov_->get_ncases());
	}

	// Set staging
	if (tp_.power) {
	  stage_ = Stage::Power;
	} else {
	  stage_ = (tp_.stage_1_permutations > 0) ? Stage::Stage1 : Stage::Stage2;
	}

	// Handle zipped input
	if (is_gzipped(tp.genotypes_path)) {
	  gt_ifs_.open(tp.genotypes_path, std::ios_base::in | std::ios_base::binary);
	  gt_streambuf.push(boost::iostreams::gzip_decompressor());
	  gt_streambuf.push(gt_ifs_);
	} else {
	  gt_ifs_.open(tp.genotypes_path, std::ios_base::in);
	  gt_streambuf.push(gt_ifs_);
	}

	std::istream gt_stream(&gt_streambuf);
	// Retrieve header line
	std::getline(gt_stream, header_);


	// Sort covariates
	cov_->sort_covariates(header_);

	permutation_ptr_ = std::make_shared<std::vector<std::vector<int32_t>>>();
	if (tp.external) { // Read in external stage 1 permutations
	  if(!check_file_exists(tp.external_path)) {
	    throw(std::runtime_error("ERROR: External permutation file doesn't exist."));
	  }
	  std::ifstream ifs(tp.external_path);
	  std::string line;
	  while(std::getline(ifs, line)) {
	    permutation_ptr_->push_back({});
	    RJBUtil::Splitter<std::string> splitter(line, " \t");
		permutation_ptr_->back().reserve(splitter.size());
		for (const auto &v : splitter) {
		  permutation_ptr_->back().emplace_back(std::stoi(v));
	    }
	  }
	  std::cerr << permutation_ptr_->size() << " external stage 1 permutations read in." << std::endl;
	  tp.stage_1_permutations = permutation_ptr_->size();
	} else {
	  // Generate permutations for stage 1
	  arma::wall_clock timer;
	  timer.tic();
	  if (tp.stage_1_permutations > 0 && !tp.alternate_permutation) {
		permute_.get_permutations(permutation_ptr_,
								  cov_->get_odds(),
								  cov_->get_ncases(),
								  tp.stage_1_permutations,
								  tp.nthreads - 1);
	  }
	  if (tp.stage_1_permutations > 0) {
		std::cerr << "Time spent generating stage 1 permutations: " << timer.toc() << std::endl;
	  } else {
		timer.toc();
	  }
	}
	// Time for 1000 permutations, 1000 samples -> 0.5s, 10000 samples-> 30s, 100000 samples 3300s

	// Permutation set output | Developer option
	if (tp.permute_set) {
	  std::ofstream pset_ofs(*tp.permute_set);
	  std::ofstream lr_ofs(*tp.permute_set + ".lr");
	  for (const auto &p : *permutation_ptr_) {
		for (const auto &v : p) {
		  pset_ofs << v << "\t";
		}
		pset_ofs << "\n";
	  }
	  pset_ofs.close();

	  cov_->get_fitted().t().print(lr_ofs);
	  lr_ofs.close();
	  // Finish if there are no stage 2 permutations.
	  if (tp.stage_2_permutations == 0) {
		std::exit(0);
	  }
	}

	if (!tp.gene_list) {
	  // Parse if no gene_list
	  all_gene_dispatcher(gt_stream);
	} else {
	  // Parse with gene_list
	  gene_list_dispatcher(gt_stream);
	}

	// Close open files
	gt_ifs_.close();

	// Wait for queue to finish processing
	tq_.join();

	// TODO: Free permutation memory
	if (permutation_ptr_.unique()) {
	  permutation_ptr_.reset();
	}

	if (tp_.gene_list) {
	  reporter->report(tq_.get_results(), tp_);
	}

	reporter->cleanup(tp_);

	cov_.reset();
  }

private:
  // Member functions
  // Dispatch
  void all_gene_dispatcher(std::istream &gt_stream) {
	std::string line;
	std::stringstream current;

	int gene_no = 1;

	while (std::getline(gt_stream, line)) {
	  RJBUtil::Splitter<std::string> split(line, "\t");
	  RJBUtil::Splitter<std::string> vsplit(split[2], "-");
	  if (vsplit.size() < 2) {
		throw (std::logic_error("Position not formatted correctly. Should be chromosome-start-end-type."));
	  }

	  // If this line is part of the same gene
	  if (split[0] == gene_) {
		add_line(current, line, split, vsplit);
	  } else {
		// Have we read a gene yet?
		if (!gene_.empty()) {
		  if (std::any_of(nvariants_.cbegin(), nvariants_.cend(), [&](const auto &v) { return v.second > 0; })) {
			// Dispatch gene
			Gene gene_data(current, cov_->get_nsamples(), nvariants_, weight_, tp_, cov_->get_original_phenotypes());

			if (tp_.range_start && tp_.range_end) {
			  if (gene_no >= *tp_.range_start && gene_no <= *tp_.range_end) {
				single_dispatch(gene_data);
			  } else if (gene_no > *tp_.range_end) {
				break; // Stop early
			  }
			} else {
			  single_dispatch(gene_data);
			}
			// Reset for next gene
			gene_no++;
			new_gene(current, line, split, vsplit);
		  } else {
              gene_no++;
              new_gene(current, line, split, vsplit);
		  }
		} else {
		  // Setup initial gene
		  new_gene(current, line, split, vsplit);
		}
	  }
	}
	Gene gene_data(current, cov_->get_nsamples(), nvariants_, weight_, tp_, cov_->get_original_phenotypes());

	if (!gene_data.is_skippable())
	  single_dispatch(gene_data);
  }
  void single_dispatch(Gene &gene) {
	using namespace std::literals::chrono_literals;
	Task_t ta(stage_,
			  gene,
			  cov_,
			  tp_,
			  *permutation_ptr_);
	// Limit adding jobs to prevent excessive memory usage
	while (tq_.size() > tp_.nthreads - 1) {
	  std::this_thread::sleep_for(0.001s);
	}
	tq_.dispatch(std::move(ta));
	ngenes_++;
  }

  // Gene list
  auto find_gene(const std::string &gene) {
	return std::find(gene_list_.cbegin(), gene_list_.cend(), gene);
  }

  void gene_list_dispatcher(std::istream &gt_stream) {
	std::string line;
	std::stringstream current;

	while (std::getline(gt_stream, line)) {
	  RJBUtil::Splitter<std::string> split(line, "\t");
	  RJBUtil::Splitter<std::string> vsplit(split[2], "-");

	  // If this line is part of the same gene
	  if (split[0] == gene_) {
		add_line(current, line, split, vsplit);
	  } else {
	    ntranscripts_ = nvariants_.size();
		// Have we read a gene yet?
		if (!gene_.empty()) {
		  if (std::any_of(nvariants_.cbegin(), nvariants_.cend(), [&](const auto &v) { return v.second > 0; })) {
			// Dispatch gene
			Gene gene_data(current, cov_->get_nsamples(), nvariants_, weight_, tp_, cov_->get_original_phenotypes());

			if (!gene_data.is_skippable())
			  multiple_dispatch(gene_data);

			auto fit = find_gene(split[0]);
			if (fit != gene_list_.cend()) {
			  // Next gene is in list
			  gene_list_.erase(fit);
			  new_gene(current, line, split, vsplit);
			} else {
			  // Check if we're done
			  if (gene_list_.empty()) {
				return;
			  }
			  // Reset to initial state to find next gene
			  reset_gene(current);
			}
		  }
		} else {
		  // Skip until we find the first gene in our list.
		  auto fit = find_gene(split[0]);
		  if (fit == gene_list_.cend()) {
			continue;
		  } else {
			gene_list_.erase(fit);
		  }
		  // Setup initial gene
		  new_gene(current, line, split, vsplit);
		}
	  }
	}
	if (!current.str().empty()) {
	  Gene gene_data(current, cov_->get_nsamples(), nvariants_, weight_, tp_, cov_->get_original_phenotypes());

	  if (!gene_data.is_skippable())
		multiple_dispatch(gene_data);
	}
  }
  void multiple_dispatch(Gene &gene) {
	using namespace std::literals::chrono_literals;
	// Ensure we have at least one variant for a submitted gene
	if (std::any_of(nvariants_.cbegin(), nvariants_.cend(), [&](const auto &v) { return v.second > 0; })) {
	  int total_s1_perm = 0;
	  int total_s2_perm = 0;
	  int total_success = 0;

	  // Single dispatch of gene list items for power analysis
	  if (tp_.power) {
		Task_t ta(stage_,
				  gene,
				  cov_,
				  tp_,
				  tp_.success_threshold,
				  tp_.stage_1_permutations,
				  tp_.stage_2_permutations,
				  *permutation_ptr_);
		while (tq_.size() > tp_.nthreads - 1) {
		  std::this_thread::sleep_for(0.001s);
		}
		tq_.dispatch(ta);
		ngenes_++;
		return;
	  }

	  for (int i = 0; i < tq_.get_nthreads(); i++) {
		if (i == tq_.get_nthreads() - 1) {
		  Task_t ta(stage_,
					gene,
					cov_,
					tp_,
					tp_.success_threshold - total_success,
					tp_.stage_1_permutations - total_s1_perm,
					tp_.stage_2_permutations - total_s2_perm,
					*permutation_ptr_);

		  // Limit adding jobs to prevent excessive memory usage
		  while (tq_.size() > tp_.nthreads - 1) {
			std::this_thread::sleep_for(0.001s);
		  }
		  tq_.dispatch(ta);
		} else {
		  Task_t ta(stage_,
					gene,
					cov_,
					tp_,
					tp_.success_threshold / static_cast<int>(tq_.get_nthreads()),
					tp_.stage_1_permutations / static_cast<int>(tq_.get_nthreads()),
					tp_.stage_2_permutations / static_cast<int>(tq_.get_nthreads()),
					*permutation_ptr_);
		  // Add current permutations
		  total_s1_perm += tp_.stage_1_permutations / tq_.get_nthreads();
		  total_s2_perm += tp_.stage_2_permutations / tq_.get_nthreads();
		  total_success += tp_.success_threshold / tq_.get_nthreads();

		  // Limit adding jobs to prevent excessive memory usage
		  while (tq_.size() > tp_.nthreads - 1) {
			std::this_thread::sleep_for(0.001s);
		  }
		  tq_.dispatch(ta);
		}
	  }
	  ngenes_++;
	}
  }

  // Input parsing
  void new_gene(std::stringstream &ss,
				std::string &line,
				RJBUtil::Splitter<std::string> &split,
				RJBUtil::Splitter<std::string> &vsplit) {
	// Reset the read buffer
	ss.str("");
	ss.clear();

	// Reset the transcript map
	nvariants_.clear();

	// Append header
	ss << header_ << "\n";
	add_line(ss, line, split, vsplit);

	// Add gene name
	gene_ = split[0];
  }
  void reset_gene(std::stringstream &ss) {
	// Reset the read buffer
	ss.str("");
	ss.clear();

	// Reset the transcript map
	nvariants_.clear();

	// Clear gene and transcript
	gene_ = "";
	transcript_ = "";
  }

  void add_line(std::stringstream &ss,
				std::string &line,
				RJBUtil::Splitter<std::string> &split,
				RJBUtil::Splitter<std::string> &vsplit) {
	// Is variant masked?
	if (!bed_.check_variant(vsplit[0], vsplit[1])) {
	  // Variant not masked
	  ss << line << "\n";
	  // Track number of variants in each transcript
	  if (nvariants_.find(split[1]) == nvariants_.end()) {
	    nvariants_[split[1]] = 1;
	  } else {
	    nvariants_[split[1]]++;
	  }
	}
  }

  // Member variables
  TaskParams tp_;
  TaskQueue2<Operation_t, Task_t, Reporter_t> tq_;
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
