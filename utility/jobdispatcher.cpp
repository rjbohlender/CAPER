//
// Created by Bohlender,Ryan James on 10/8/18.
//

#include <chrono>

#include "jobdispatcher.hpp"

#include "main_support.hpp"
#include "reporter.hpp"
#include "filesystem.hpp"

using namespace std::chrono_literals;

JobDispatcher::JobDispatcher(TaskParams &tp, std::shared_ptr<Reporter> reporter)
	: tp_(tp), tq_(tp_.nthreads - 1, reporter, tp.verbose),
	  cov_(std::make_shared<Covariates>(tp.covariates_path, tp.ped_path, tp.linear)) {
  // Initialize bed and weights
  if (tp.gene_list)
	gene_list_ = RJBUtil::Splitter<std::string>(*tp.gene_list, ",");
  if (tp.bed)
	bed_ = Bed(*tp.bed);
  if (tp.weight)
	weight_ = Weight(*tp.weight);

  // Set staging
  stage_ = (tp_.stage_1_permutations > 0) ? Stage::Stage1 : Stage::Stage2;

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
  std::cerr << "Time spent generating stage 1 permutations: " << timer.toc() << std::endl;
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

  if (tp.gene_list)
	Reporter(tq_.get_results(), tp_);

  cov_.reset();
}

void JobDispatcher::all_gene_dispatcher(std::istream &gt_stream) {
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

		  if(tp_.range_start && tp_.range_end) {
		    if(gene_no >= *tp_.range_start && gene_no <= *tp_.range_end) {
			  single_dispatch(gene_data);
		    } else if(gene_no > *tp_.range_end) {
		      break; // Stop early
		    }
		  } else {
			single_dispatch(gene_data);
		  }
		  // Reset for next gene
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

void JobDispatcher::single_dispatch(Gene &gene) {
  TaskArgs ta(stage_,
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

auto JobDispatcher::find_gene(const std::string &gene) {
  return std::find(gene_list_.cbegin(), gene_list_.cend(), gene);
}

void JobDispatcher::gene_list_dispatcher(std::istream &gt_stream) {
  std::string line;
  std::stringstream current;

  while (std::getline(gt_stream, line)) {
	RJBUtil::Splitter<std::string> split(line, "\t");
	RJBUtil::Splitter<std::string> vsplit(split[2], "-");

	// If this line is part of the same gene
	if (split[0] == gene_) {
	  add_line(current, line, split, vsplit);
	} else {
	  // Have we read a gene yet?
	  if (!gene_.empty()) {
		if (std::any_of(nvariants_.cbegin(), nvariants_.cend(), [&](const auto &v) { return v.second > 0; })) {
		  // Dispatch gene
		  Gene gene_data(current, cov_->get_nsamples(), nvariants_, weight_, tp_);

		  if(!gene_data.is_skippable())
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

	if(!gene_data.is_skippable())
	  multiple_dispatch(gene_data);
  }
}

void JobDispatcher::multiple_dispatch(Gene &gene) {
  // Ensure we have at least one variant for a submitted gene
  if (std::any_of(nvariants_.cbegin(), nvariants_.cend(), [&](const auto &v) { return v.second > 0; })) {
	int total_s1_perm = 0;
	int total_s2_perm = 0;
	int total_success = 0;

	for (int i = 0; i < tq_.get_nthreads(); i++) {
	  if (i == tq_.get_nthreads() - 1) {
		TaskArgs ta(stage_,
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
		TaskArgs ta(stage_,
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

void JobDispatcher::new_gene(std::stringstream &ss,
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

  // New transcript
  ntranscripts_++;

  // Add gene name
  gene_ = split[0];
}

void JobDispatcher::reset_gene(std::stringstream &ss) {
  // Reset the read buffer
  ss.str("");
  ss.clear();

  // Reset the transcript map
  nvariants_.clear();

  // Clear gene and transcript
  gene_ = "";
  transcript_ = "";
}

void JobDispatcher::add_line(std::stringstream &ss,
							 std::string &line,
							 RJBUtil::Splitter<std::string> &split,
							 RJBUtil::Splitter<std::string> &vsplit) {
  // Is variant masked?
  if (!bed_.check_variant(vsplit[0], vsplit[1])) {
	// Variant not masked
	ss << line << "\n";
	// Track number of variants in each transcript
	if (split[1] == transcript_) {
	  nvariants_[transcript_]++;
	} else {
	  transcript_ = split[1];
	  nvariants_[transcript_] = 1;
	  ntranscripts_++;
	}
  }
}
