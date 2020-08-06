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
#include "../data/permutation.hpp"
#include "../data/bed.hpp"
#include "../data/weight.hpp"
#include "taskqueue.hpp"
#include "../carva/carvaop.hpp"

#include "reporter.hpp"
#include "filesystem.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

template<typename Operation_t, typename Task_t, typename Reporter_t>
class JobDispatcher {
public:
  JobDispatcher(TaskParams &tp, std::shared_ptr<Reporter_t> reporter)
	  : tp_(tp), tq_(tp.nthreads - 1, reporter, tp),
		cov_(std::make_shared<Covariates>(tp.covariates_path, tp.ped_path, tp.linear)) {
    if (tp_.seed) {
      permute_ = Permute(*tp_.seed);
    }
	// Initialize bed and weights
	if (tp_.gene_list) {
	  gene_list_ = RJBUtil::Splitter<std::string>(*tp_.gene_list, ",");
	}
	if (tp_.bed) {
	  bed_ = Bed(*tp_.bed);
	}
	if (tp_.weight) {
	  weight_ = Weight(*tp_.weight);
	}

	// Update case/control count for reporter
	if(!tp_.power) {
	  reporter->set_ncases(cov_->get_ncases());
	  reporter->set_ncontrols(cov_->get_nsamples() - cov_->get_ncases());
	}

	// Set staging
	if (tp_.power) {
	  stage_ = Stage::Power;
	} else {
	  stage_ = Stage::Stage1;
	}

	// Handle zipped input
	if (is_gzipped(tp_.genotypes_path)) {
	  gt_ifs_.open(tp_.genotypes_path, std::ios_base::in | std::ios_base::binary);
	  gt_streambuf.push(boost::iostreams::gzip_decompressor());
	  gt_streambuf.push(gt_ifs_);
	} else {
	  gt_ifs_.open(tp_.genotypes_path, std::ios_base::in);
	  gt_streambuf.push(gt_ifs_);
	}

	std::istream gt_stream(&gt_streambuf);
	// Retrieve header line
	std::getline(gt_stream, header_);

	// Cleanup previous permutations if looping because we append to the file
	if(tp_.permute_set) {
	  if(check_file_exists(*tp_.permute_set)) {
	    std::remove((*tp_.permute_set).c_str());
	  }
	}

	// Sort covariates
	cov_->sort_covariates(header_);

	permutation_ptr_ = std::make_shared<std::vector<std::vector<int8_t>>>();
	if(!tp_.max_perms) {
	  if (tp_.external) { // Read in external stage 1 permutations
		if (!check_file_exists(tp_.external_path)) {
		  throw (std::runtime_error("ERROR: External permutation file doesn't exist."));
		}
		std::ifstream ifs(tp_.external_path);
		std::string line;
		while (std::getline(ifs, line)) {
		  permutation_ptr_->push_back({});
		  permutation_ptr_->back().reserve(line.size());
		  for (const auto &v : line) {
			if (v == '0' || v == '1') {
			  permutation_ptr_->back().emplace_back(v - '0');
			}
		  }
		}
		std::cerr << permutation_ptr_->size() << " external stage 1 permutations read in." << std::endl;
		tp_.nperm = permutation_ptr_->size();
	  } else if(tp.nocovadj) {
		// Generate permutations for stage 1
		arma::wall_clock timer;
		timer.tic();
		if (tp_.nperm > 0) {
		  *permutation_ptr_ = std::vector<std::vector<int8_t>>(tp_.nperm, std::vector<int8_t>(cov_->get_nsamples(), 0));
		  for(auto &p : *permutation_ptr_) {
			for(int i = 0; i < cov_->get_ncases(); i++) { // Assign cases
			  p[i] = 1;
			}
			permute_.fisher_yates(p, permute_.sto);
		  }
		}
		if (tp_.nperm > 0 && tp_.verbose) {
		  std::cerr << "Time spent generating unadjusted permutations: " << timer.toc() << std::endl;
		}
	  } else {
		// Generate permutations for stage 1
		arma::wall_clock timer;
		timer.tic();
		if (tp_.nperm > 0) {
		  permute_.generate_permutations(permutation_ptr_,
										 cov_->get_odds(),
										 cov_->get_ncases(),
										 tp_.nperm,
										 tp_.nthreads - 1,
										 tp_.bin_epsilon);
		}
		if (tp_.nperm > 0 && tp_.verbose) {
		  std::cerr << "Time spent generating adjusted permutations: " << timer.toc() << std::endl;
		}
	  }
	  // Time for 1000 permutations, 1000 samples -> 0.5s, 10000 samples-> 30s, 100000 samples 3300s

	  // Permutation set output | Developer option
	  if (tp_.permute_set) {
		std::ofstream pset_ofs(*tp_.permute_set);
		for (const auto &p : *permutation_ptr_) {
		  for (const auto &v : p) {
			pset_ofs << static_cast<char>(v + '0');
		  }
		  pset_ofs << "\n";
		}
		pset_ofs.close();
	  }

	  if (!tp_.gene_list) {
		// Parse if no gene_list
		all_gene_dispatcher(gt_stream);
	  } else {
		// Parse with gene_list
		gene_list_dispatcher(gt_stream);
	  }

	  gt_ifs_.close();
	} else {
	  arma::uword remaining = *tp_.max_perms;
	  while(remaining > 0) {
		if (tp_.external) {
		  if (!check_file_exists(tp_.external_path)) {
			throw (std::runtime_error("ERROR: External permutation file doesn't exist."));
		  }
		  std::ifstream ifs(tp_.external_path);
		  ifs.seekg(external_pos); // Doesn't move on first loop.
		  std::string line;
		  unsigned lineno = 0;
		  while (std::getline(ifs, line)) {
			permutation_ptr_->push_back({});
			permutation_ptr_->back().reserve(line.size());
			for (const auto &v : line) {
			  if (v == '0' || v == '1') {
				permutation_ptr_->back().emplace_back(v - '0');
			  }
			}
			lineno++;
			if(lineno >= tp_.nperm) {
			  break;
			}
		  }
		  external_pos = ifs.tellg(); // Record our position and continue form there
		  remaining -= lineno;
		} else if(tp.nocovadj) {
		  // Generate permutations for stage 1
		  arma::wall_clock timer;
		  timer.tic();
		  if (tp_.nperm > 0) {
		    *permutation_ptr_ = std::vector<std::vector<int8_t>>(std::min(tp_.nperm, remaining), std::vector<int8_t>(cov_->get_nsamples(), 0));
		    for(auto &p : *permutation_ptr_) {
		      for(int i = 0; i < cov_->get_ncases(); i++) { // Assign cases
		        p[i] = 1;
		      }
		      permute_.fisher_yates(p, permute_.sto);
		    }
		  }
		  if (tp_.nperm > 0 && tp_.verbose) {
			std::cerr << "Time spent generating unadjusted permutations: " << timer.toc() << std::endl;
		  }
		  remaining -= std::min(tp_.nperm, remaining);
		} else {
		  // Generate permutations for stage 1
		  arma::wall_clock timer;
		  timer.tic();
		  if (tp_.nperm > 0) {
			permute_.generate_permutations(permutation_ptr_,
										   cov_->get_odds(),
										   cov_->get_ncases(),
										   std::min(tp_.nperm, remaining),
										   tp_.nthreads - 1,
										   tp_.bin_epsilon);
		  }
		  if (tp_.nperm > 0 && tp_.verbose) {
			std::cerr << "Time spent generating adjusted permutations: " << timer.toc() << std::endl;
		  }
		  remaining -= std::min(tp_.nperm, remaining);
		}
		if (tp_.permute_set) {
		  std::ofstream pset_ofs(*tp_.permute_set, std::ios::app | std::ios::out);
		  for (const auto &p : *permutation_ptr_) {
			for (const auto &v : p) {
			  pset_ofs << static_cast<char>(v + '0');
			}
			pset_ofs << "\n";
		  }
		  pset_ofs.close();
		}
		if(gt_ifs_.is_open()) {
		  if (!tp_.gene_list) {
			// Parse if no gene_list
			all_gene_dispatcher(gt_stream);
		  } else {
			// Parse with gene_list
			gene_list_dispatcher(gt_stream);
		  }
		  gt_ifs_.close();
		  tq_.wait(); // Need to wait here until all jobs are done, otherwise permutations will update during processing
		  if(tq_.continue_.size() == 0) { // Terminate if all jobs are done.
		    remaining = 0;
		  }
		} else { // Successive loops will not need to process the gene stream
		  tq_.redispatch(); // Re-dispatch each job. No need to update anything because the permutations are kept in a pointer.
		  auto cur = std::chrono::system_clock::now();
		  std::time_t cur_time = std::chrono::system_clock::to_time_t(cur);
		  std::cerr << "Remaining permutations: " << remaining << " at " << std::ctime(&cur_time);
		  tq_.wait(); // Need to wait here until all jobs are done, otherwise permutations will update during processing
		  if(tq_.continue_.size() == 0) { // Terminate if all jobs are done.
			remaining = 0;
		  }
		}
	  }
	}

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
			Gene gene_data(current, cov_->get_nsamples(), nvariants_, weight_, tp_);

			if(!gene_data.is_skippable()) {
              if (tp_.range_start && tp_.range_end) {
                if (gene_no >= *tp_.range_start && gene_no <= *tp_.range_end) {
                  single_dispatch(gene_data);
                } else if (gene_no > *tp_.range_end) {
                  break; // Stop early
                }
              } else {
                single_dispatch(gene_data);
              }
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
	Gene gene_data(current, cov_->get_nsamples(), nvariants_, weight_, tp_);

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
			Gene gene_data(current, cov_->get_nsamples(), nvariants_, weight_, tp_);

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
	  Gene gene_data(current, cov_->get_nsamples(), nvariants_, weight_, tp_);

	  if (!gene_data.is_skippable())
		multiple_dispatch(gene_data);
	}
  }
  void multiple_dispatch(Gene &gene) {
	using namespace std::literals::chrono_literals;
	// Ensure we have at least one variant for a submitted gene
	if (std::any_of(nvariants_.cbegin(), nvariants_.cend(), [&](const auto &v) { return v.second > 0; })) {
	  long total_perm = 0;
	  long total_success = 0;
	  long perm_step = tp_.nperm / (tp_.nthreads - 1);
	  long succ_step = tp_.success_threshold / (tp_.nthreads - 1);
	  long max_loops = 1;
	  if(tp_.max_perms) {
	    max_loops = *tp_.max_perms / tp_.nperm;
	  }

	  // Single dispatch of gene list items for power analysis
	  if (tp_.power) {
		Task_t ta(stage_,
				  gene,
				  cov_,
				  tp_,
				  tp_.success_threshold,
				  tp_.nperm,
				  0,
				  tp_.nperm,
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
					tp_.nperm - total_perm,
					i * perm_step,
					max_loops * (tp_.nperm - total_perm),
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
					succ_step,
					perm_step,
					i * perm_step,
					max_loops * perm_step,
					*permutation_ptr_);
		  // Add current permutations
		  total_perm += perm_step;
		  total_success += succ_step;

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
	if (vsplit[2] != vsplit[1]) {
      if (!bed_.check_variant(vsplit[0], std::make_pair(vsplit[1], vsplit[2]))) {
        // Variant not masked
        ss << line << "\n";
        // Track number of variants in each transcript
        if (nvariants_.find(split[1]) == nvariants_.end()) {
          nvariants_[split[1]] = 1;
        } else {
          nvariants_[split[1]]++;
        }
      }
	} else {
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
  }

  // Member variables
  TaskParams tp_;
  TaskQueue<Operation_t, Task_t, Reporter_t> tq_;
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
  size_t external_pos = 0;

  std::shared_ptr<Covariates> cov_;
  std::shared_ptr<std::vector<std::vector<int8_t>>> permutation_ptr_;
};

#endif //PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP
