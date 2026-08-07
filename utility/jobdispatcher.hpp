//
// Created by Bohlender,Ryan James on 10/8/18.
//

#ifndef PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP
#define PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP

#include <algorithm>
#include <chrono>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "../caper/caperop.hpp"
#include "../caper/capertask.hpp"
#include "../data/bed.hpp"
#include "../data/covariates.hpp"
#include "../data/permutation.hpp"
#include "../data/weight.hpp"
#include "reporter.hpp"
#include "taskqueue.hpp"

#include "../data/filter.hpp"
#include "filesystem.hpp"
#include "filevalidator.hpp"
#include "reporter.hpp"

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

template <typename Operation_t, typename Task_t, typename Reporter_t>
class JobDispatcher {
public:
  JobDispatcher(TaskParams &tp, std::shared_ptr<Reporter_t> reporter)
      : tp_(tp), tq_(tp.nthreads - 1, reporter, tp),
        cov_(std::make_shared<Covariates>(tp)) {
    if (tp_.seed) {
      permute_ = Permute(*tp_.seed);
    }
    // Initialize bed and weights
    if (tp_.gene_list) {
      RJBUtil::Splitter<std::string> gene_split(*tp_.gene_list, ",");
      for (const auto &gene_token : gene_split) {
        gene_list_.emplace(RJBUtil::Splitter<std::string>::str(gene_token));
      }
    }
    if (tp_.input_paths.empty()) {
      tp_.input_paths.push_back(tp_.input_path);
    }
    if (tp_.input_paths.size() > 1 && gene_list_.size() != 1) {
      throw std::runtime_error(
          "Multiple input files require exactly one gene specified with --genes.");
    }
    if (tp_.bed) {
      bed_ = Bed(*tp_.bed);
    }
    if (tp_.weight) {
      weight_ = Weight(*tp_.weight);
    }
    Filter filter(tp_.whitelist_path);

    // Update case/control count for reporter
    if (!tp_.power) {
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
    boost::iostreams::file_source gt_file_{tp_.input_path};
    if (is_gzipped(tp_.input_path)) {
      gt_ifs_.push(boost::iostreams::gzip_decompressor());
      gt_ifs_.push(gt_file_);
    } else if(is_zstd(tp_.input_path)) {
      gt_ifs_.push(boost::iostreams::zstd_decompressor());
      gt_ifs_.push(gt_file_);
    } else {
      gt_ifs_.push(gt_file_);
    }

    // Retrieve header line
    if (!std::getline(gt_ifs_, header_)) {
      throw std::runtime_error("Unable to read matrix header from: " +
                               tp_.input_path);
    }

    if (tp_.input_paths.size() > 1) {
      prepare_multi_input_gene(*gene_list_.begin());
    }

    // Cleanup previous permutations if looping because we append to the file
    if (tp_.permute_set) {
      if (check_file_exists(*tp_.permute_set)) {
        std::remove((*tp_.permute_set).c_str());
      }
    }

    // Sort covariates
    cov_->sort_covariates(header_);

    permutation_ptr_ = std::make_shared<std::vector<std::vector<int8_t>>>();
    if (!tp_.max_perms) {
      if (tp_.external) { // Read in external stage 1 permutations
        if (!check_file_exists(tp_.external_path)) {
          throw(std::runtime_error(
              "ERROR: External permutation file doesn't exist."));
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
        std::cerr << permutation_ptr_->size()
                  << " external stage 1 permutations read in." << std::endl;
        tp_.nperm = permutation_ptr_->size();
      } else if (tp.nocovadj) {
        // Generate permutations for stage 1
        arma::wall_clock timer;
        timer.tic();
        if (tp_.nperm > 0) {
          *permutation_ptr_ = std::vector<std::vector<int8_t>>(
              tp_.nperm, std::vector<int8_t>(cov_->get_nsamples(), 0));
          for (auto &p : *permutation_ptr_) {
            for (int i = 0; i < cov_->get_ncases(); i++) { // Assign cases
              p[i] = 1;
            }
            Permute::fisher_yates(p, permute_.sto);
          }
        }
        if (tp_.nperm > 0 && tp_.verbose) {
          std::cerr << "Time spent generating unadjusted permutations: "
                    << timer.toc() << std::endl;
        }
      } else {
        // Generate permutations for stage 1
        arma::wall_clock timer;
        timer.tic();
        if (tp_.nperm > 0) {
          permute_.generate_permutations(permutation_ptr_, cov_->get_odds(),
                                         cov_->get_ncases(), tp_.nperm,
                                         tp_.nthreads - 1, tp_.bin_epsilon);
        }
        if (tp_.nperm > 0 && tp_.verbose) {
          std::cerr << "Time spent generating adjusted permutations: "
                    << timer.toc() << std::endl;
        }
      }
      // Time for 1000 permutations, 1000 samples -> 0.5s, 10000 samples-> 30s,
      // 100000 samples 3300s

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

      dispatch_input(filter);

    } else {
      arma::uword remaining = *tp_.max_perms;
      while (remaining > 0) {
        if (tp_.external) {
          if (!check_file_exists(tp_.external_path)) {
            throw(std::runtime_error(
                "ERROR: External permutation file doesn't exist."));
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
            if (lineno >= tp_.nperm) {
              break;
            }
          }
          external_pos =
              ifs.tellg(); // Record our position and continue from there
          remaining -= lineno;
        } else if (tp.nocovadj) {
          // Generate permutations for stage 1
          arma::wall_clock timer;
          timer.tic();
          if (tp_.nperm > 0) {
            *permutation_ptr_ = std::vector<std::vector<int8_t>>(
                std::min(tp_.nperm, remaining),
                std::vector<int8_t>(cov_->get_nsamples(), 0));
            for (auto &p : *permutation_ptr_) {
              for (int i = 0; i < cov_->get_ncases(); i++) { // Assign cases
                p[i] = 1;
              }
              Permute::fisher_yates(p, permute_.sto);
            }
          }
          if (tp_.nperm > 0 && tp_.verbose) {
            std::cerr << "Time spent generating unadjusted permutations: "
                      << timer.toc() << std::endl;
          }
          remaining -= std::min(tp_.nperm, remaining);
        } else {
          // Generate permutations for stage 1
          arma::wall_clock timer;
          timer.tic();
          if (tp_.nperm > 0) {
            permute_.generate_permutations(permutation_ptr_, cov_->get_odds(),
                                           cov_->get_ncases(),
                                           std::min(tp_.nperm, remaining),
                                           tp_.nthreads - 1, tp_.bin_epsilon);
          }
          if (tp_.nperm > 0 && tp_.verbose) {
            std::cerr << "Time spent generating adjusted permutations: "
                      << timer.toc() << std::endl;
          }
          remaining -= std::min(tp_.nperm, remaining);
        }
        if (tp_.permute_set) {
          std::ofstream pset_ofs(*tp_.permute_set,
                                 std::ios::app | std::ios::out);
          for (const auto &p : *permutation_ptr_) {
            for (const auto &v : p) {
              pset_ofs << static_cast<char>(v + '0');
            }
            pset_ofs << "\n";
          }
          pset_ofs.close();
        }
#if 1
        // First loop
        if (remaining >= *tp_.max_perms - tp_.nperm) {
	        dispatch_input(filter);
	        tq_.wait(); // Need to wait here until all jobs are done, otherwise
	        // permutations will update during processing
	        if (tq_.continue_.size() == 0) { // Terminate if all jobs are done.
		        remaining = 0;
	        }
        } else {
          // Subsequent loops
	        tq_.redispatch(); // Re-dispatch each job. No need to update anything
	        // because the permutations are kept in a pointer.
	        auto cur = std::chrono::system_clock::now();
	        std::time_t cur_time = std::chrono::system_clock::to_time_t(cur);
	        std::cerr << "Remaining permutations: " << remaining << " at "
	                  << std::ctime(&cur_time);
	        tq_.wait(); // Need to wait here until all jobs are done, otherwise
	        // permutations will update during processing
	        if (tq_.continue_.size() == 0) { // Terminate if all jobs are done.
		        remaining = 0;
	        }
        }
#else
        if (gt_ifs_.good()) {
          if (!tp_.gene_list) {
            // Parse if no gene_list
            all_gene_dispatcher(gt_ifs_, filter);
          } else {
            // Parse with gene_list
            gene_list_dispatcher(gt_ifs_, filter);
          }
          tq_.wait(); // Need to wait here until all jobs are done, otherwise
                      // permutations will update during processing
          if (tq_.continue_.size() == 0) { // Terminate if all jobs are done.
            remaining = 0;
          }
        } else { // Successive loops will not need to process the gene stream
          tq_.redispatch(); // Re-dispatch each job. No need to update anything
                            // because the permutations are kept in a pointer.
          auto cur = std::chrono::system_clock::now();
          std::time_t cur_time = std::chrono::system_clock::to_time_t(cur);
          std::cerr << "Remaining permutations: " << remaining << " at "
                    << std::ctime(&cur_time);
          tq_.wait(); // Need to wait here until all jobs are done, otherwise
                      // permutations will update during processing
          if (tq_.continue_.size() == 0) { // Terminate if all jobs are done.
            remaining = 0;
          }
        }
#endif
      }
    }

    // Wait for queue to finish processing
    tq_.join();

    // TODO: Free permutation memory
    if (!tp_.external && permutation_ptr_.use_count() == 1) {
      permutation_ptr_.reset();
    }

    if (tp_.gene_list) {
      reporter->report(tq_.get_results(), tp_);
    }

    if (tp_.method == "VAAST" && !tp_.gene_list && !tp_.power) {
      reporter->vaast_sample_index_map(cov_->get_samples());
    }
    reporter->cleanup(tp_);

    cov_.reset();
  }

private:
  struct MatrixHeader {
    std::vector<std::string> metadata_columns;
    std::vector<std::string> samples;
  };

  struct CollectedLine {
    std::string transcript;
    std::string variant_key;
    std::string line;
    std::string path;
    int line_number;
  };

  // Member functions
  // Dispatch
  void dispatch_input(Filter &filter) {
    if (tp_.input_paths.size() > 1) {
      multi_file_single_gene_dispatcher(filter);
    } else if (tp_.gene_list) {
      gene_list_dispatcher(gt_ifs_, filter);
    } else {
      all_gene_dispatcher(gt_ifs_, filter);
    }
  }

  static MatrixHeader parse_matrix_header(const std::string &header,
                                          const std::string &path) {
    RJBUtil::Splitter<std::string> fields(header, "\t", 0, true);
    const auto first_sample = static_cast<size_t>(Indices::first);
    if (fields.size() <= first_sample) {
      throw std::runtime_error("Matrix header contains no sample columns: " +
                               path);
    }

    MatrixHeader result;
    result.metadata_columns.reserve(first_sample);
    result.samples.reserve(fields.size() - first_sample);
    for (size_t i = 0; i < first_sample; ++i) {
      result.metadata_columns.push_back(fields.str(i));
    }

    std::unordered_set<std::string> seen;
    for (size_t i = first_sample; i < fields.size(); ++i) {
      auto sample = fields.str(i);
      if (sample.empty()) {
        throw std::runtime_error(
            "Matrix header contains an empty sample ID: " + path);
      }
      if (!seen.insert(sample).second) {
        throw std::runtime_error("Matrix header contains duplicate sample ID '" +
                                 sample + "': " + path);
      }
      result.samples.push_back(std::move(sample));
    }
    return result;
  }

  static std::vector<size_t>
  validate_header_and_build_mapping(const MatrixHeader &canonical,
                                    const MatrixHeader &current,
                                    const std::string &path) {
    if (current.metadata_columns != canonical.metadata_columns) {
      throw std::runtime_error(
          "Matrix metadata columns differ from the first input file: " +
          path);
    }

    std::unordered_set<std::string> canonical_samples(
        canonical.samples.begin(), canonical.samples.end());
    std::unordered_set<std::string> current_samples(current.samples.begin(),
                                                    current.samples.end());
    std::vector<std::string> missing;
    std::vector<std::string> unexpected;

    for (const auto &sample : canonical.samples) {
      if (!current_samples.contains(sample)) {
        missing.push_back(sample);
      }
    }
    for (const auto &sample : current.samples) {
      if (!canonical_samples.contains(sample)) {
        unexpected.push_back(sample);
      }
    }

    if (!missing.empty() || !unexpected.empty() ||
        current.samples.size() != canonical.samples.size()) {
      std::ostringstream message;
      message << "Matrix samples differ from the first input file: " << path;
      if (!missing.empty()) {
        message << ". Missing samples: ";
        for (size_t i = 0; i < missing.size(); ++i) {
          if (i > 0) message << ", ";
          message << missing[i];
        }
      }
      if (!unexpected.empty()) {
        message << ". Unexpected samples: ";
        for (size_t i = 0; i < unexpected.size(); ++i) {
          if (i > 0) message << ", ";
          message << unexpected[i];
        }
      }
      throw std::runtime_error(message.str());
    }

    std::unordered_map<std::string, size_t> current_positions;
    current_positions.reserve(current.samples.size());
    for (size_t i = 0; i < current.samples.size(); ++i) {
      current_positions.emplace(current.samples[i], i);
    }

    std::vector<size_t> mapping;
    mapping.reserve(canonical.samples.size());
    for (const auto &sample : canonical.samples) {
      mapping.push_back(current_positions.at(sample));
    }
    return mapping;
  }

  static std::string reorder_genotypes(
      std::string_view genotypes, const std::vector<size_t> &mapping,
      const std::string &path, int line_number) {
    if (genotypes.size() != mapping.size()) {
      throw std::runtime_error(
          "Genotype count differs from the matrix header in " + path +
          " at line " + std::to_string(line_number));
    }

    std::string reordered(mapping.size(), '\0');
    for (size_t canonical_index = 0; canonical_index < mapping.size();
         ++canonical_index) {
      reordered[canonical_index] = genotypes[mapping[canonical_index]];
    }
    return reordered;
  }

  static std::string rebuild_matrix_line(
      RJBUtil::Splitter<std::string> &fields,
      const std::string &genotypes) {
    std::ostringstream result;
    const auto first_sample = static_cast<size_t>(Indices::first);
    for (size_t i = 0; i < first_sample; ++i) {
      if (i > 0) result << '\t';
      result << fields[i];
    }
    result << '\t' << genotypes;
    return result.str();
  }

  void prepare_multi_input_gene(const std::string &target_gene) {
    const auto canonical =
        parse_matrix_header(header_, tp_.input_paths.front());
    std::unordered_map<std::string, std::pair<std::string, int>> variants;

    for (const auto &path : tp_.input_paths) {
      boost::iostreams::filtering_istream input;
      boost::iostreams::file_source file(path);
      if (!file.is_open()) {
        throw std::runtime_error("Unable to open genotype matrix: " + path);
      }
      if (is_gzipped(path)) {
        input.push(boost::iostreams::gzip_decompressor());
      } else if (is_zstd(path)) {
        input.push(boost::iostreams::zstd_decompressor());
      }
      input.push(file);

      std::string file_header;
      if (!std::getline(input, file_header)) {
        throw std::runtime_error("Unable to read matrix header from: " + path);
      }
      const auto current = parse_matrix_header(file_header, path);
      const auto mapping =
          validate_header_and_build_mapping(canonical, current, path);
      FileValidator validator;
      validator.set_matrix_header(file_header);

      std::string line;
      int line_number = 1;
      while (std::getline(input, line)) {
        ++line_number;
        RJBUtil::Splitter<std::string> prefix(line, "\t", 7, true);
        if (prefix.size() <= static_cast<size_t>(Indices::gene)) {
          throw std::runtime_error("Truncated matrix line in " + path +
                                   " at line " +
                                   std::to_string(line_number));
        }
        if (prefix.str(static_cast<int>(Indices::gene)) != target_gene) {
          continue;
        }

        RJBUtil::Splitter<std::string> fields(line, "\t", 11, true);
        validator.validate_matrix_line(fields, line_number);
        const auto reordered = reorder_genotypes(
            fields[static_cast<int>(Indices::first)], mapping, path,
            line_number);

        std::ostringstream key;
        key << fields[static_cast<int>(Indices::chrom)] << '\x1f'
            << fields[static_cast<int>(Indices::start)] << '\x1f'
            << fields[static_cast<int>(Indices::end)] << '\x1f'
            << fields[static_cast<int>(Indices::ref)] << '\x1f'
            << fields[static_cast<int>(Indices::alt)] << '\x1f'
            << fields[static_cast<int>(Indices::transcript)];
        const auto variant_key = key.str();
        const auto [found, inserted] =
            variants.emplace(variant_key, std::make_pair(path, line_number));
        if (!inserted) {
          throw std::runtime_error(
              "Duplicate variant/transcript entry for gene " + target_gene +
              " in " + path + " at line " + std::to_string(line_number) +
              "; first seen in " + found->second.first + " at line " +
              std::to_string(found->second.second));
        }

        multi_input_lines_.push_back(
            {fields.str(static_cast<int>(Indices::transcript)), variant_key,
             rebuild_matrix_line(fields, reordered), path, line_number});
      }
    }

    if (multi_input_lines_.empty()) {
      throw std::runtime_error("Gene " + target_gene +
                               " was not found in any input matrix.");
    }
    std::stable_sort(multi_input_lines_.begin(), multi_input_lines_.end(),
                     [](const CollectedLine &lhs, const CollectedLine &rhs) {
                       return lhs.transcript < rhs.transcript;
                     });
    multi_input_gene_ = target_gene;
  }

  void multi_file_single_gene_dispatcher(Filter &filter) {
    std::stringstream current;
    reset_gene(current);
    current << header_ << '\n';
    for (auto &collected : multi_input_lines_) {
      RJBUtil::Splitter<std::string> fields(collected.line, "\t", 11, true);
      add_line(current, collected.line, fields);
    }
    gene_ = multi_input_gene_;

    if (!std::any_of(nvariants_.cbegin(), nvariants_.cend(),
                     [](const auto &value) { return value.second > 0; })) {
      return;
    }
    Gene gene_data(current, cov_, cov_->get_nsamples(), nvariants_, weight_,
                   tp_, filter);
    if (!gene_data.is_skippable()) {
      multiple_dispatch(gene_data);
    }
  }

  void all_gene_dispatcher(std::istream &gt_stream, Filter &filter) {
    std::string line;
    std::stringstream current;
    std::unordered_set<std::string> previous_genes;
    int lineno = -1;
    FileValidator fv;
    fv.set_matrix_header(header_);

    int gene_no = 1;

    while (std::getline(gt_stream, line)) {
      lineno++;
      RJBUtil::Splitter<std::string> split(line, "\t", 11, true);
      fv.validate_matrix_line(split, lineno);

      // If this line is part of the same gene
      if (split[static_cast<int>(Indices::gene)] == gene_) {
        add_line(current, line, split);
      } else {
        auto gene_token = split.str(static_cast<int>(Indices::gene));
        if (previous_genes.contains(gene_token)) {
          throw std::runtime_error(
              "Gene list must be sorted by gene name. Gene " + gene_token +
              " appears again on line " + std::to_string(lineno) +
              " of the gene stream. Please sort the gene stream by gene name "
              "and transcript.");
        } else {
          previous_genes.insert(gene_token);
        }
        // Have we read a gene yet?
        if (!gene_.empty()) {
          if (std::any_of(nvariants_.cbegin(), nvariants_.cend(),
                          [&](const auto &v) { return v.second > 0; })) {
            // Dispatch gene
            Gene gene_data(current, cov_, cov_->get_nsamples(), nvariants_,
                           weight_, tp_, filter);

            if (!gene_data.is_skippable()) {
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
            new_gene(current, line, split);
          } else {
            gene_no++;
            new_gene(current, line, split);
          }
        } else {
          // Setup initial gene
          new_gene(current, line, split);
        }
      }
    }
    Gene gene_data(current, cov_, cov_->get_nsamples(), nvariants_, weight_,
                   tp_, filter);

    if (!gene_data.is_skippable())
      single_dispatch(gene_data);
  }

  void single_dispatch(Gene &gene) {
    Task_t ta(stage_, gene, cov_, tp_, *permutation_ptr_);
    tq_.wait_for_space(tp_.nthreads - 1);
    tq_.dispatch(std::move(ta));
    ngenes_++;
  }

  void gene_list_dispatcher(std::istream &gt_stream, Filter &filter) {
    std::string line;
    std::stringstream current;
    int lineno = -1;
    FileValidator fv;
    fv.set_matrix_header(header_);

    while (std::getline(gt_stream, line)) {
      lineno++;
      RJBUtil::Splitter<std::string> split(line, "\t", 7);

      // If this line is part of the same gene
      if (split[static_cast<int>(Indices::gene)] == gene_) {
        // Split the rest of the way
        split = RJBUtil::Splitter<std::string>(line, "\t");
        fv.validate_matrix_line(split, lineno);
        add_line(current, line, split);
      } else {
        // Have we read a gene yet?
        if (!gene_.empty()) {
          if (std::any_of(nvariants_.cbegin(), nvariants_.cend(),
                          [&](const auto &v) { return v.second > 0; })) {
            // Dispatch gene
            Gene gene_data(current, cov_, cov_->get_nsamples(), nvariants_,
                           weight_, tp_, filter);

            if (!gene_data.is_skippable())
              multiple_dispatch(gene_data);

            auto gene_name = split.str(static_cast<int>(Indices::gene));
            if (gene_list_.contains(gene_name)) {
              // Next gene is in list
              split = RJBUtil::Splitter<std::string>(line, "\t");
              fv.validate_matrix_line(split, lineno);
              gene_list_.erase(gene_name);
              new_gene(current, line, split);
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
          auto gene_name = split.str(static_cast<int>(Indices::gene));
          if (!gene_list_.contains(gene_name)) {
            continue;
          }
          gene_list_.erase(gene_name);
          split = RJBUtil::Splitter<std::string>(line, "\t");
          fv.validate_matrix_line(split, lineno);
          // Setup initial gene
          new_gene(current, line, split);
        }
      }
    }
    if (!current.str().empty()) {
      Gene gene_data(current, cov_, cov_->get_nsamples(), nvariants_, weight_,
                     tp_, filter);

      if (!gene_data.is_skippable())
        multiple_dispatch(gene_data);
    }
  }

  void multiple_dispatch(Gene &gene) {
    // Ensure we have at least one variant for a submitted gene
    if (std::any_of(nvariants_.cbegin(), nvariants_.cend(),
                    [&](const auto &v) { return v.second > 0; })) {

      long n_workers = tq_.get_nthreads();
      if (n_workers <= 1) {
          Task_t ta(stage_, gene, cov_, tp_, tp_.success_threshold, tp_.nperm, 0,
                    tp_.nperm, *permutation_ptr_);
          tq_.wait_for_space(tp_.nthreads - 1);
          tq_.dispatch(std::move(ta));
          ngenes_++;
          return;
      }

      long total_perm = 0;
      long total_success = 0;
      long perm_step = tp_.nperm / n_workers;
      long succ_step = tp_.success_threshold / n_workers;
      long max_loops = 1;
      if (tp_.max_perms) {
        max_loops = std::ceil(*tp_.max_perms / (double)tp_.nperm);
      }


      // Single dispatch of gene list items for power analysis
      if (tp_.power) {
        Task_t ta(stage_, gene, cov_, tp_, tp_.success_threshold, tp_.nperm, 0,
                  tp_.nperm, *permutation_ptr_);
        tq_.wait_for_space(tp_.nthreads - 1);
        tq_.dispatch(std::move(ta));
        ngenes_++;
        return;
      }

      for (int i = 0; i < n_workers; i++) {
        if (i == n_workers - 1) {
          Task_t ta(stage_,
                    gene,
                    cov_,
                    tp_,
                    tp_.success_threshold - total_success,
                    tp_.nperm - total_perm,
                    i * perm_step,
                    max_loops * (tp_.nperm - total_perm),
                    *permutation_ptr_);

          tq_.wait_for_space(n_workers);
          tq_.dispatch(std::move(ta));
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

          tq_.wait_for_space(n_workers);
          tq_.dispatch(std::move(ta));
        }
      }
      ngenes_++;
    }
  }

  // Input parsing
  void new_gene(std::stringstream &ss, std::string &line,
                RJBUtil::Splitter<std::string> &split) {
    // Reset the read buffer
    ss.str("");
    ss.clear();

    // Reset the transcript map
    nvariants_.clear();

    // Append header
    ss << header_ << "\n";
    add_line(ss, line, split);

    // Add gene name
    gene_ = split[static_cast<int>(Indices::gene)];
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

  void add_line(std::stringstream &ss, std::string &line,
                RJBUtil::Splitter<std::string> &split) {
    // Build variant id once and check mask status
    const auto chrom = split[static_cast<int>(Indices::chrom)];
    const auto start = split[static_cast<int>(Indices::start)];
    const auto end = split[static_cast<int>(Indices::end)];
    const auto ref = split[static_cast<int>(Indices::ref)];
    const auto alt = split[static_cast<int>(Indices::alt)];

    std::string variant_id;
    variant_id.reserve(chrom.size() + start.size() + end.size() + ref.size() +
                       alt.size() + 4);
    variant_id.append(chrom.data(), chrom.size());
    variant_id.push_back(',');
    variant_id.append(start.data(), start.size());
    variant_id.push_back(',');
    variant_id.append(end.data(), end.size());
    variant_id.push_back(',');
    variant_id.append(ref.data(), ref.size());
    variant_id.push_back(',');
    variant_id.append(alt.data(), alt.size());

    if (!bed_.check_variant(variant_id)) {
      // Variant not masked
      ss << line << "\n";
      // Track number of variants in each transcript
      auto transcript = split.str(static_cast<int>(Indices::transcript));
      ++nvariants_[transcript];
    }
  }

  // Member variables
  TaskParams tp_;
  std::unordered_set<std::string> gene_list_;
  TaskQueue<Operation_t, Task_t, Reporter_t> tq_;
  boost::iostreams::filtering_istream gt_ifs_;
  Bed bed_;
  Weight weight_;
  Permute permute_;

  // Gene parsing
  std::string header_;
  std::string gene_;
  std::string transcript_;
  std::string multi_input_gene_;
  std::vector<CollectedLine> multi_input_lines_;
  Stage stage_;
  std::unordered_map<std::string, arma::uword> nvariants_;

  // Counters
  arma::uword ngenes_ = 0;
  size_t external_pos = 0;

  std::shared_ptr<Covariates> cov_;
  std::shared_ptr<std::vector<std::vector<int8_t>>> permutation_ptr_;
};

#endif // PERMUTE_ASSOCIATE_JOBDISPATCHER_HPP
