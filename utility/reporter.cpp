//
// Created by Bohlender,Ryan James on 10/11/18.
//

#include "reporter.hpp"
#include "filesystem.hpp"
#include "math.hpp"

#include <boost/format.hpp>

const std::set<std::string> Reporter::pvalue_methods_{"SKAT", "SKATO", "CMC",
                                                      "RVT1", "RVT2"};

Reporter::Reporter(TaskParams &tp)
    : ncases_(0), ncontrols_(0), method_(tp.method), pvalues_(tp.analytic),
      gene_list_(tp.gene_list), print_testable_(tp.testable) {
  std::string header;

  if (!check_directory_exists(tp.output_path)) {
    throw(std::runtime_error("Output path is invalid."));
  }

  // If using VAAST, output VAAST file
  if (tp.method == "VAAST") {
    vaast_path_ss << tp.output_path << "/" << tp.method;
    if (tp.group_size > 0)
      vaast_path_ss << boost::format(".g%1$d") % tp.group_size;
    if (tp.testable)
      vaast_path_ss << ".check_testability";
    if (tp.biallelic)
      vaast_path_ss << ".biallelic";
    if (tp.range_start && tp.range_end) {
      vaast_path_ss << "." << *tp.range_start;
      vaast_path_ss << "." << *tp.range_end;
    }
    vaast_path_ss << ".vaast";
    vaast_file_ = std::ofstream(vaast_path_ss.str());
  }

  simple_path_ss << tp.output_path << "/" << tp.method;
  if (tp.group_size > 0 && tp.method == "VAAST")
    simple_path_ss << boost::format(".g%1$d") % tp.group_size;
  if (tp.testable)
    simple_path_ss << ".check_testability";
  if (tp.biallelic && tp.method == "VAAST")
    simple_path_ss << ".biallelic";
  if (tp.range_start && tp.range_end) {
    simple_path_ss << "." << *tp.range_start;
    simple_path_ss << "." << *tp.range_end;
  }
  simple_path_ss << ".simple";
  simple_path_tmp_ss << simple_path_ss.str() << ".tmp";

  simple_file_tmp_ = std::ofstream(simple_path_tmp_ss.str());
  if (!simple_file_tmp_.good()) {
    throw(std::runtime_error("Simple file failed to open for writing.\n"));
  }
  simple_file_tmp_ << std::setw(25) << std::left << "Gene"
                   << " ";
  simple_file_tmp_ << std::setw(20) << "Transcript";
  simple_file_tmp_ << std::setw(20) << "Test_Statistic";
  simple_file_tmp_ << std::setw(20) << "Exact_P";
  simple_file_tmp_ << std::setw(20) << "Empirical_P";
  simple_file_tmp_ << std::setw(20) << "Empirical_P_CI";
  simple_file_tmp_ << std::setw(20) << "Empirical_MidP";
  simple_file_tmp_ << std::setw(20) << "Empirical_MidP_CI";
  simple_file_tmp_ << std::setw(20) << "MGIT_P";
  simple_file_tmp_ << std::setw(20) << "MGIT_MIDP";
  simple_file_tmp_ << std::setw(20) << "Successes";
  simple_file_tmp_ << std::setw(20) << "MGIT_Successes";
  simple_file_tmp_ << std::setw(20) << "MGIT_MIDP_Successes";
  simple_file_tmp_ << std::setw(20) << "Permutations" << std::endl;

  detail_path_ss << tp.output_path << "/" << tp.method;
  if (tp.group_size > 0 && tp.method == "VAAST")
    detail_path_ss << boost::format(".g%1$d") % tp.group_size;
  if (tp.testable)
    detail_path_ss << ".check_testability";
  if (tp.biallelic && tp.method == "VAAST")
    detail_path_ss << ".biallelic";
  if (tp.range_start && tp.range_end) {
    detail_path_ss << "." << *tp.range_start;
    detail_path_ss << "." << *tp.range_end;
  }
  detail_path_ss << ".detail";

  detail_file_ = std::ofstream(detail_path_ss.str());
  if (!detail_file_.good()) {
    throw(std::runtime_error("Detail file failed to open for writing.\n"));
  }

  if (tp.qtl) {
    header = "Gene\tTranscripts\tVariant\tScore\tAF";
  } else {
    header = "Gene\tTranscripts\tVariant\tScore\tWeight\tAF\tcase_ref\tcase_"
             "alt\tcontrol_ref\tcontrol_alt\tcase_list\tcontrol_list";
  }
  detail_file_ << header << std::endl;
}

Reporter::Reporter(std::vector<CAPERTask> &res, TaskParams &tp)
    : ncases_(0), ncontrols_(0), method_(tp.method), pvalues_(tp.analytic),
      gene_list_(tp.gene_list), print_testable_(tp.testable) {
  if (!check_directory_exists(tp.output_path)) {
    throw(std::runtime_error("Output path is invalid."));
  }

  // Normal execution
  simple_path_ss << tp.output_path << "/" << tp.method;
  if (tp.group_size > 0 && tp.method == "VAAST")
    simple_path_ss << boost::format(".g%1$d") % tp.group_size;
  if (tp.testable)
    simple_path_ss << ".check_testability";
  if (tp.biallelic && tp.method == "VAAST")
    simple_path_ss << ".biallelic";
  if (tp.range_start && tp.range_end) {
    simple_path_ss << "." << *tp.range_start;
    simple_path_ss << "." << *tp.range_end;
  }
  simple_path_ss << ".simple";

  detail_path_ss << tp.output_path << "/" << tp.method;
  if (tp.group_size > 0 && tp.method == "VAAST")
    detail_path_ss << boost::format(".g%1$d") % tp.group_size;
  if (tp.testable)
    detail_path_ss << ".check_testability";
  if (tp.biallelic && tp.method == "VAAST")
    detail_path_ss << ".biallelic";
  if (tp.range_start && tp.range_end) {
    detail_path_ss << "." << *tp.range_start;
    detail_path_ss << "." << *tp.range_end;
  }
  detail_path_ss << ".detail";

  simple_file_tmp_ = std::ofstream(simple_path_ss.str());
  detail_file_ = std::ofstream(detail_path_ss.str());

  if (!simple_file_tmp_.good()) {
    throw(std::runtime_error("Simple file failed to open for writing.\n"));
  }
  if (!detail_file_.good()) {
    throw(std::runtime_error("Detail file failed to open for writing.\n"));
  }
  // Extract results
  extract_results(res, tp);

  // Write output
  if (tp.gene_list) {
    report_simple(tp);
    if (!tp.no_detail) {
      report_detail(res, tp);
    }
  }
}

auto Reporter::report(std::vector<CAPERTask> &&res, TaskParams &tp) -> void {
  // Extract results
  extract_results(res, tp);

  // Write output
  if (tp.gene_list) {
    report_simple(tp);
    if (!tp.no_detail) {
      report_detail(res, tp);
      if (tp.method == "VAAST") {
        report_vaast(res, tp);
      }
    }
  }
}

auto Reporter::cleanup(TaskParams &tp) -> void { sort_simple(tp); }

auto Reporter::extract_results(std::vector<CAPERTask> &tq_results,
                               TaskParams &tp) -> void {
  if (gene_list_) {
    // Combine results
    for (auto &task : tq_results) {
      for (auto &result : task.results) {
        if (results_.find(result.second.gene) != results_.end()) {
          if (results_[result.second.gene].find(result.second.transcript) !=
              results_[result.second.gene].end()) {
            results_[result.second.gene][result.second.transcript].combine(
                result.second, tp);
          } else {
            results_[result.second.gene][result.second.transcript] =
                result.second;
            vaast_.emplace(std::make_pair(
                result.second.transcript,
                task.gene.get_vaast()[result.second.transcript]));
          }
        } else {
          results_[result.second.gene][result.second.transcript] =
              result.second;
          details_.push_back(task.gene.get_detail());
          vaast_.emplace(
              std::make_pair(result.second.transcript,
                             task.gene.get_vaast()[result.second.transcript]));
        }
      }
    }
    // Recalculate MGIT
    recalculate_mgit(results_);
  } else {
    for (auto &task : tq_results) {
      for (auto &result : task.results) {
        results_[result.second.gene][result.second.transcript] = result.second;
        details_.push_back(task.gene.get_detail());
      }
    }
  }

#if 0
  // Sort extracted results
  if(pvalue_methods_.find(method_) != pvalue_methods_.end() && (method_ != "SKAT" || tp.nperm == 0 )) {
    // Sort by midp and asymptotic pvalue
    std::sort(results_.begin(), results_.end(), [](auto &a, auto &b) {
      if(a.empirical_midp != b.empirical_midp) {
        return a.empirical_midp < b.empirical_midp;
      } else {
        return a.original < b.original;
      }
    });
  } else {
    // Sort by midp and score
    std::sort(results_.begin(), results_.end(), [](auto &a, auto &b) {
      if(a.empirical_midp != b.empirical_midp) {
        return a.empirical_midp < b.empirical_midp;
      } else {
        return a.original > b.original;
      }
    });
  }
#endif
}

auto Reporter::recalculate_mgit(
    std::map<std::string, std::map<std::string, Result>> &results) -> void {
  for (auto &g : results) {
    // Skip MGIT if only a single transcript
    if (g.second.size() == 1) {
      for (auto &v : g.second) {
        v.second.mgit_p = v.second.empirical_p;
        v.second.mgit_midp = v.second.empirical_midp;
        v.second.mgit_successes = v.second.successes;
        v.second.mgit_midp_successes = v.second.mid_successes;
      }
      continue;
    }

    // Shorthand
    unsigned long n = g.second.size();
    arma::uword max_perm = 0;
    arma::uword i;
    double successes = 0;
    double midp_successes = 0;

    // Get max_perm
    for (const auto &tr : g.second) {
      if (tr.second.permutations > max_perm) {
        max_perm = tr.second.permutations;
      }
    }

    arma::mat mgit_pval_mat = arma::mat(max_perm, n);
    i = 0;
    for (auto &tr : g.second) {
      int m = tr.second.permuted.size();

      arma::vec permuted = arma::conv_to<arma::vec>::from(tr.second.permuted);
      arma::vec pvals;

      if (pvalues_) {
        pvals = rank(permuted, "ascend");
      } else {
        pvals = rank(permuted, "descend");
      }

      pvals /= permuted.n_rows;

      try {
        mgit_pval_mat.col(i) = pvals;
      } catch (const std::logic_error &e) {
        std::cerr << "n_row: " << mgit_pval_mat.n_rows
                  << " n_col: " << mgit_pval_mat.n_cols << "\n";
        throw(e);
      }
      i++;
    }
    arma::vec mgit_pval_dist_ = arma::min(mgit_pval_mat, 1);

    double min_p = std::numeric_limits<double>::max();
    for (auto &tr : g.second) {
      if (tr.second.empirical_p < min_p) {
        min_p = tr.second.empirical_p;
      }
    }

    unsigned long m = mgit_pval_dist_.n_rows; // Total permutations
    for (auto &tr : g.second) {
      successes = arma::find(mgit_pval_dist_ <= min_p).eval().n_rows;
      midp_successes = arma::find(mgit_pval_dist_ < min_p).eval().n_rows;
      midp_successes += arma::find(mgit_pval_dist_ == min_p).eval().n_rows / 2.;

      // Store multi-transcript p-value
      tr.second.mgit_p = (1.0 + successes) / (1.0 + m);
      tr.second.mgit_successes = static_cast<int>(successes);
      tr.second.mgit_midp = (1.0 + midp_successes) / (1.0 + m);
      tr.second.mgit_midp_successes = static_cast<int>(midp_successes);
    }
  }
}

auto Reporter::recalculate_mgit(
    std::unordered_map<std::string, Result> &results) -> void {
  // Skip MGIT if only a single transcript
  if (results.size() == 1) {
    for (auto &v : results) {
      v.second.mgit_p = v.second.empirical_p;
      v.second.mgit_midp = v.second.empirical_midp;
      v.second.mgit_successes = v.second.successes;
      v.second.mgit_midp_successes = v.second.mid_successes;
    }
    return;
  }

  // Shorthand
  unsigned long n = results.size(); // Number of transcripts
  arma::uword max_perm = 0;
  arma::uword i, j, k;
  double successes = 0;
  double midp_successes = 0;

  // Get max_perm
  for (const auto &tr : results) {
    if (tr.second.permutations > max_perm) {
      max_perm = tr.second.permutations;
    }
  }

  arma::mat mgit_pval_mat = arma::mat(max_perm, n);
  i = 0;
  for (auto &tr : results) {
    int m = tr.second.permuted.size();
    if (m == 0) {
      continue;
    }

    arma::vec permuted = arma::conv_to<arma::vec>::from(tr.second.permuted);
    arma::vec pvals;

    // SKATO and SKAT return pvalues so reverse success criteria
    if (pvalues_) {
      pvals = rank(permuted, "ascend");
    } else {
      pvals = rank(permuted, "descend");
    }

    pvals /= permuted.n_rows;

    try {
      mgit_pval_mat.col(i) = pvals;
    } catch (const std::logic_error &e) {
      std::cerr << "n_row: " << mgit_pval_mat.n_rows
                << " n_col: " << mgit_pval_mat.n_cols << "\n";
      throw(e);
    }
    i++;
  }
  arma::vec mgit_pval_dist_ = arma::min(mgit_pval_mat, 1);

  double min_p = std::numeric_limits<double>::max();
  for (auto &tr : results) {
    if (tr.second.empirical_p < min_p) {
      min_p = tr.second.empirical_p;
    }
  }

  for (auto &tr : results) {
    unsigned long m = mgit_pval_dist_.n_rows; // Total permutations

    successes = arma::find(mgit_pval_dist_ <= min_p).eval().n_rows;
    midp_successes = arma::find(mgit_pval_dist_ < min_p).eval().n_rows;
    midp_successes += arma::find(mgit_pval_dist_ == min_p).eval().n_rows / 2.;

    // Store multi-transcript p-value
    tr.second.mgit_p = (1.0 + successes) / (1.0 + m);
    tr.second.mgit_successes = static_cast<int>(successes);
    tr.second.mgit_midp = (1.0 + midp_successes) / (1.0 + m);
    tr.second.mgit_midp_successes = static_cast<int>(midp_successes);
  }
}

auto Reporter::report_simple(TaskParams &tp) -> void {
  // Print header and formatted results
  double permutation_mean = 0;
  double permutation_variance = 0;

  for (auto &resmap : results_) {
    for (auto &result : resmap.second) {
      if (result.second.successes < tp.success_threshold) {
        if (std::find(unfinished_.begin(), unfinished_.end(),
                      result.second.gene) == unfinished_.end())
          unfinished_.push_back(result.second.gene);
      }
      permutation_mean += result.second.permutations;
    }
  }
  permutation_mean /= results_.size();
  for (auto &resmap : results_) {
    for (auto &result : resmap.second) {
      permutation_variance +=
          std::pow(result.second.permutations - permutation_mean, 2) /
          results_.size();
    }
  }
  // Report results
  int rank = 1;
  for (auto &resmap : results_) {
    for (auto &result : resmap.second) {
      result.second.rank = rank;
      write_to_stream(simple_file_tmp_, result.second);
      rank++;
    }
  }

  simple_file_tmp_.flush();
  simple_file_tmp_.close();

  std::cerr << "Permutation mean: " << permutation_mean << std::endl;
  std::cerr << "Permutation sd: " << std::sqrt(permutation_variance)
            << std::endl;
  std::cerr << "Transcripts submitted: " << results_.size() << std::endl;
}

auto Reporter::report_detail(std::vector<CAPERTask> &res, TaskParams &tp)
    -> void {
  int i = 0; // For each gene
  for (auto &v : details_) {
    detail_file_ << v;
    // Print sample / index map at the end
    if (i == res.size() - 1) {
      detail_file_ << "## Sample Index Map" << std::endl;
      int j = 0;
      for (auto &s : res[0].gene.get_samples()) {
        detail_file_ << "#\t" << s << "\t" << j << std::endl;
        j++;
      }
    }
    i++;
  }
}

auto Reporter::report_vaast(std::vector<CAPERTask> &res, TaskParams &tp)
    -> void {
  // Header information
  vaast_file_ << "## PA_VERSION\t0.0" << std::endl;
  vaast_file_ << "## COMMAND\t" << tp.full_command << std::endl;

  for (auto &resmap : results_) {
    for (auto &result : resmap.second) {
      vaast_file_ << vaast_[result.second.transcript];

      vaast_file_ << "SCORE: "
                  << boost::format("%1$.2f") % result.second.original
                  << std::endl;
      vaast_file_ << "genome_permutation_p: " << result.second.empirical_p
                  << std::endl;
      vaast_file_ << "genome_permutation_p_ci: "
                  << result.second.empirical_ci.first << ","
                  << result.second.empirical_ci.second << std::endl;
      vaast_file_ << "num_permutations: " << result.second.permutations
                  << std::endl;
      vaast_file_ << "total_success: " << result.second.successes << std::endl;
    }
  }
  vaast_sample_index_map(res);
}
void Reporter::vaast_sample_index_map(const std::vector<CAPERTask> &res) {
  // Print sample / index map at the end
  vaast_file_ << "## Sample Index Map" << std::endl;
  int j = 0;
  if (results_.empty()) {
    return;
  }
  for (const auto &s : res[0].gene.get_samples()) {
    vaast_file_ << "#\t" << s << "\t" << j << std::endl;
    j++;
  }
}

void Reporter::vaast_sample_index_map(std::vector<std::string> &&samples) {
  // Print sample / index map at the end
  vaast_file_ << "## Sample Index Map" << std::endl;
  int j = 0;
  for (const auto &s : samples) {
    vaast_file_ << "#\t" << s << "\t" << j << std::endl;
    j++;
  }
}

auto Reporter::write_to_stream(std::ostream &os, Result &res) -> void {
  if (print_testable_) {
    if (res.testable) {
      os << res;
    }
  } else {
    os << res;
  }
}

auto Reporter::sync_write_simple(
    std::unordered_map<std::string, Result> &results, const TaskParams &tp)
    -> void {
  std::lock_guard<std::mutex> lock(lock_); // Acquire lock

  recalculate_mgit(results);

  for (const auto &v : results) {
    if (!v.second.skippable) {
      if (v.second.successes < tp.success_threshold) {
        if (std::find(unfinished_.begin(), unfinished_.end(), v.second.gene) ==
            unfinished_.end())
          unfinished_.push_back(v.second.gene);
      }
    }
  }
  if (tp.top_only) {
    // Find the most significant result.
    std::unique_ptr<Result> topres;
    for (auto &r : results) {
      if (r.second.skippable) {
        continue;
      }
      if (topres == nullptr) {
        topres = std::make_unique<Result>(r.second);
      } else {
        if (topres->empirical_p > r.second.empirical_p) {
          topres = std::make_unique<Result>(r.second);
        }
      }
    }
    if (topres != nullptr) {
      if (print_testable_) {
        if (topres->testable) {
          simple_file_tmp_ << *topres;
        }
      } else {
        simple_file_tmp_ << *topres;
      }
    }
    return;
  }

  if (print_testable_) {
    for (auto &tr : results) {
      if (tr.second.testable && !tr.second.skippable) {
        simple_file_tmp_ << tr.second;
      }
    }
  } else {
    for (auto &tr : results) {
      if (tr.second.skippable) {
        continue;
      }
      simple_file_tmp_ << tr.second;
    }
  }

  simple_file_tmp_.flush();
}

auto Reporter::sync_write_detail(const std::string &d, bool gene_testable)
    -> void {
  std::lock_guard<std::mutex> lock(lock_); // Acquire lock

  if (print_testable_) {
    if (gene_testable) {
      detail_file_ << d;
    }
  } else {
    detail_file_ << d;
  }

  detail_file_.flush();
}

auto Reporter::sort_simple(const TaskParams &tp) -> void {
  std::lock_guard<std::mutex> lock(lock_);
  simple_file_ = std::ofstream(simple_path_ss.str());
  if (!simple_file_.good()) {
    throw(std::runtime_error("Simple file failed to open for writing.\n"));
  }
  simple_file_tmp_.close();
  std::ifstream ifs(simple_path_tmp_ss.str());
  std::string line;

  unsigned long rank = 1;

  simple_file_ << "# Command: " << tp.full_command << std::endl;
  simple_file_ << "# Cases: " << ncases_ << std::endl;
  simple_file_ << "# Controls: " << ncontrols_ << std::endl;
  simple_file_ << std::setw(20) << std::left << "Rank";
  simple_file_ << std::setw(25) << "Gene"
               << " ";
  simple_file_ << std::setw(20) << "Transcript" << " ";
  simple_file_ << std::setw(30) << "Test_Statistic";
  // simple_file_ << std::setw(20) << "Exact_P";
  simple_file_ << std::setw(20) << "Empirical_P";
  simple_file_ << std::setw(20) << "Empirical_P_CI";
  simple_file_ << std::setw(20) << "Empirical_MidP";
  simple_file_ << std::setw(20) << "Empirical_MidP_CI";
  simple_file_ << std::setw(20) << "MGIT_P";
  simple_file_ << std::setw(20) << "MGIT_MIDP";
  simple_file_ << std::setw(20) << "Successes";
  simple_file_ << std::setw(20) << "MGIT_Successes";
  simple_file_ << std::setw(20) << "MGIT_MIDP_Successes";
  simple_file_ << std::setw(20) << "Permutations" << std::endl;

  std::vector<ResultLine> results;
  unsigned long lineno = 0;

  while (std::getline(ifs, line)) {
    if (lineno == 0) {
      lineno++;
      continue;
    }
    RJBUtil::Splitter<std::string> splitter(line, " \t");
    RJBUtil::Splitter<std::string> emp_ci_splitter(splitter[4], ",");
    RJBUtil::Splitter<std::string> emp_midci_splitter(splitter[6], ",");

    try {
      ResultLine rs =
          ResultLine{splitter[0], splitter[1], std::stod(splitter[2]),
                     // std::stod(splitter[3]),
                     std::stod(splitter[3]),
                     std::make_pair(std::stod(emp_ci_splitter[0]),
                                    std::stod(emp_ci_splitter[1])),
                     std::stod(splitter[5]),
                     std::make_pair(std::stod(emp_midci_splitter[0]),
                                    std::stod(emp_midci_splitter[1])),
                     std::stod(splitter[7]), std::stod(splitter[8]),
                     std::stoul(splitter[9]), std::stoul(splitter[10]),
                     std::stoul(splitter[11]), std::stoul(splitter[12])};

      for (int i = 12; i < splitter.size(); i++) {
        rs.stats.push_back(splitter[i]);
      }
      results.push_back(rs);
    } catch (std::exception &e) {
      std::cerr << "Error parsing line: " << line << std::endl
                << "On line number: " << lineno << std::endl;
      throw(e);
    }

    lineno++;
  }

  ifs.close();

  // Sort results and write them out
  if (!tp.analytic) {
    std::sort(results.begin(), results.end(), [](ResultLine &a, ResultLine &b) {
      return a.empirical_p < b.empirical_p;
    });
  } else {
    std::sort(results.begin(), results.end(), [](ResultLine &a, ResultLine &b) {
      return a.original < b.original;
    });
  }

  for (const auto &rs : results) {
    std::stringstream ci;
    std::stringstream midci;

    ci << rs.empirical_p_ci.first << "," << rs.empirical_p_ci.second;
    midci << rs.empirical_midp_ci.first << "," << rs.empirical_midp_ci.second;

    simple_file_ << std::setw(20) << std::left << rank << " ";
    simple_file_ << std::setw(25) << rs.gene << " ";
    simple_file_ << std::setw(20) << rs.transcript << " ";
    simple_file_ << std::setw(30) << std::setprecision(15) << rs.original;
    // simple_file_ << std::setw(20) << rs.exact_p;
    simple_file_ << std::setw(20) << std::setprecision(8) << rs.empirical_p;
    simple_file_ << std::setw(20) << ci.str();
    simple_file_ << std::setw(20) << rs.empirical_midp;
    simple_file_ << std::setw(20) << midci.str();
    simple_file_ << std::setw(20) << rs.mgit;
    simple_file_ << std::setw(20) << rs.mgit_midp;
    simple_file_ << std::setw(20) << rs.successes;
    simple_file_ << std::setw(20) << rs.mgit_successes;
    simple_file_ << std::setw(20) << rs.mgit_midp_successes;
    simple_file_ << std::setw(20) << rs.permutations;
    for (const auto &v : rs.stats) {
      simple_file_ << std::setw(30) << v;
    }
    simple_file_ << std::endl;

    rank++;
  }

  simple_file_.flush();
  simple_file_.close();

  // Delete tmp file
  std::remove(simple_path_tmp_ss.str().c_str());

  if (!unfinished_.empty() && tp.nperm > 0) {
    // Print command to run unfinished
    std::stringstream uf_ss;
    uf_ss << tp.program_path << " ";
    uf_ss << "-i " << tp.input_path << " ";
    uf_ss << "-c " << tp.covariates_path << " ";
    uf_ss << "-o " << tp.output_path << " ";
    uf_ss << "-p " << tp.ped_path << " ";
    if (tp.bed) {
      uf_ss << "-b " << *tp.bed << " ";
    }
    if (tp.weight) {
      uf_ss << "-w " << *tp.weight << " ";
    }
    uf_ss << "--nperm " << tp.nperm << " ";
    uf_ss << "--max_perms "
          << (tp.max_perms ? *tp.max_perms * 10 : tp.nperm * 10) << " ";
    uf_ss << "-m " << tp.method << " ";
    uf_ss << "-t " << tp.nthreads << " ";
    if (!tp.verbose) {
      uf_ss << "-q ";
    }
    if (tp.a != 1 || tp.b != 25) {
      uf_ss << "--beta_weights " << tp.a << "," << tp.b << " ";
    }
    if (tp.group_size > 0) {
      uf_ss << " -g " << tp.group_size << " ";
    }
    if (tp.maf < 0.05) {
      uf_ss << "-r " << tp.maf << " ";
    }
    if (tp.mac < std::numeric_limits<unsigned long long>::max()) {
      uf_ss << "--mac " << tp.mac << " ";
    }
    if (tp.no_detail) {
      uf_ss << "--nodetail ";
    }
    if (tp.alternate_grouping) {
      uf_ss << "--alternate_grouping ";
    }
    if (tp.soft_maf_filter != 0.5) {
      uf_ss << "--soft_maf_filter " << tp.soft_maf_filter << " ";
    }

    uf_ss << "-l ";
    for (int i = 0; i < unfinished_.size(); i++) {
      if (i == unfinished_.size() - 1) {
        uf_ss << unfinished_[i];
      } else {
        uf_ss << unfinished_[i] << ",";
      }
    }
    std::cerr << "Some genes did not reach the success threshold. Run the "
                 "following command to check those genes."
              << std::endl;
    std::cerr << uf_ss.str() << std::endl;
  }
}

auto Reporter::set_ncases(int ncases) -> void { ncases_ = ncases; }

auto Reporter::set_ncontrols(int ncontrols) -> void { ncontrols_ = ncontrols; }

auto Reporter::sync_write_vaast(CAPERTask &ct, const TaskParams &tp) -> void {
  std::lock_guard<std::mutex> lock(lock_); // Acquire lock
  if (tp.method != "VAAST") {
    return;
  }

  for (auto &ts : ct.gene.get_vaast()) {
    double z2 = std::pow(1.96, 2);
    double n = ct.results[ts.first].permutations;
    double ns = ct.results[ts.first].successes;
    double nf =
        ct.results[ts.first].permutations - ct.results[ts.first].successes;
    double p = (ns + 1.) / (n + 1.);
    double sp = (p + z2 / (2 * (n + 1.))) / (1. + z2 / (n + 1.));
    double ci = 1.96 / (1. + z2 / n) *
                std::sqrt((p * (1 - p) / (n + 1.) + z2 / (4 * n * n)));
    vaast_file_ << ts.second;
    vaast_file_ << "SCORE: "
                << boost::format("%1$.2f") % ct.results[ts.first].original
                << std::endl;
    vaast_file_ << "genome_permutation_p: " << p << std::endl;
    vaast_file_ << "genome_permutation_p_ci: " << ((sp - ci > 0) ? sp - ci : 0)
                << "," << ((sp + ci > 1) ? 1 : sp + ci) << std::endl;
    vaast_file_ << "num_permutations: " << ct.results[ts.first].permutations
                << std::endl;
    vaast_file_ << "total_success: " << ct.results[ts.first].successes
                << std::endl;
  }

  vaast_file_.flush();
}

PowerReporter::PowerReporter(TaskParams &tp) {
  std::stringstream power_path_ss;

  if (tp.method == "VAAST" && tp.group_size > 0) {
    power_path_ss << tp.output_path << "/" << tp.method << ".g" << tp.group_size
                  << ".power";
  } else {
    power_path_ss << tp.output_path << "/" << tp.method << ".power";
  }
  power_file_ = std::ofstream(power_path_ss.str());

  // Initial header write
  power_file_ << std::setw(15) << "Gene"
              << " ";
  power_file_ << std::setw(15) << "Transcript";
  power_file_ << std::setw(15) << "Method";
  power_file_ << std::setw(15) << "Ncases";
  power_file_ << std::setw(15) << "Ncontrols";
  power_file_ << std::setw(15) << "Successes";
  power_file_ << std::setw(15) << "Bootstraps";
  power_file_ << std::setw(15) << "Ratio";
  power_file_ << std::setw(15) << "Alpha";
  power_file_ << std::endl;
}

PowerReporter::PowerReporter(std::vector<PowerTask> &res, TaskParams &tp) {
  std::stringstream power_path_ss;

  if (tp.method == "VAAST" && tp.group_size > 0) {
    power_path_ss << tp.output_path << "/" << tp.method << ".g" << tp.group_size
                  << ".power";
  } else {
    power_path_ss << tp.output_path << "/" << tp.method << ".power";
  }
  power_file_ = std::ofstream(power_path_ss.str());

  // Initial header write
  power_file_ << std::setw(15) << "Gene"
              << " ";
  power_file_ << std::setw(15) << "Transcript";
  power_file_ << std::setw(15) << "Method";
  power_file_ << std::setw(15) << "Ncases";
  power_file_ << std::setw(15) << "Ncontrols";
  power_file_ << std::setw(15) << "Successes";
  power_file_ << std::setw(15) << "Bootstraps";
  power_file_ << std::setw(15) << "Ratio";
  power_file_ << std::setw(15) << "Alpha";
  power_file_ << std::endl;

  if (tp.gene_list) {
    report_power(res, tp);
  }
}

auto PowerReporter::report(std::vector<PowerTask> &resv, TaskParams &tp)
    -> void {
  report_power(resv, tp);
}

auto PowerReporter::report_power(std::vector<PowerTask> &resv, TaskParams &tp)
    -> void {
  for (const auto &pt : resv) {
    for (const auto &pr : pt.power_res_) {
      power_file_ << std::setw(25) << pr.gene << " ";
      power_file_ << std::setw(20) << pr.transcript;
      power_file_ << std::setw(20) << pr.method;
      power_file_ << std::setw(20) << pr.successes;
      power_file_ << std::setw(20) << pr.bootstraps;
      power_file_ << std::setw(20) << pr.ratio;
      power_file_ << std::setw(20) << pr.alpha;
      power_file_ << std::endl;
    }
  }
}

auto PowerReporter::sync_write_power(std::vector<PowerRes> &prv) -> void {
  std::unique_lock<std::mutex> lock(lock_); // Acquire lock

  for (const auto &pr : prv) {
    power_file_ << std::setw(25) << pr.gene << " ";
    power_file_ << std::setw(20) << pr.transcript;
    power_file_ << std::setw(20) << pr.method;
    power_file_ << std::setw(20) << pr.successes;
    power_file_ << std::setw(20) << pr.bootstraps;
    power_file_ << std::setw(20) << pr.ratio;
    power_file_ << std::setw(20) << pr.alpha;
    power_file_ << std::endl;
  }

  lock.unlock();
}

auto PowerReporter::set_ncases(int ncases) -> void { ncases_ = ncases; }

auto PowerReporter::set_ncontrols(int ncontrols) -> void {
  ncontrols_ = ncontrols;
}

auto PowerReporter::cleanup(TaskParams &tp) -> void { return; }

CAESEReporter::CAESEReporter(TaskParams &tp) {
  std::stringstream caese_path_ss;

  if (tp.method == "VAAST" && tp.group_size > 0) {
    caese_path_ss << tp.output_path << "/" << tp.method << ".g" << tp.group_size
                  << ".caese";
  } else {
    caese_path_ss << tp.output_path << "/" << tp.method << ".caese";
  }
  caese_file_ = std::ofstream(caese_path_ss.str());

  // Initial header write
  caese_file_ << std::setw(25) << std::left << "Gene"
              << " ";
  caese_file_ << std::setw(20) << std::left << "Transcript";
  caese_file_ << std::setw(20) << std::left << "OR_estimate";
  caese_file_ << std::setw(20) << std::left << "OR_CI_low";
  caese_file_ << std::setw(20) << std::left << "OR_CI_high";
  caese_file_ << std::setw(20) << std::left << "P_value";
  caese_file_ << std::setw(25) << std::left << "case_carrier";
  caese_file_ << std::setw(25) << std::left << "case_non_carrier";
  caese_file_ << std::setw(25) << std::left << "control_carrier";
  caese_file_ << std::setw(25) << std::left << "control_non_carrier";
  caese_file_ << std::endl;
}

CAESEReporter::CAESEReporter(std::vector<CAESETask> &res, TaskParams &tp) {
  std::stringstream caese_path_ss;

  if (tp.method == "VAAST" && tp.group_size > 0) {
    caese_path_ss << tp.output_path << "/" << tp.method << ".g" << tp.group_size
                  << ".caese";
  } else {
    caese_path_ss << tp.output_path << "/" << tp.method << ".caese";
  }
  caese_file_ = std::ofstream(caese_path_ss.str());

  // Initial header write
  caese_file_ << std::setw(25) << std::left << "Gene"
              << " ";
  caese_file_ << std::setw(20) << std::left << "Transcript";
  caese_file_ << std::setw(20) << std::left << "OR_estimate";
  caese_file_ << std::setw(20) << std::left << "OR_CI_low";
  caese_file_ << std::setw(20) << std::left << "OR_CI_high";
  caese_file_ << std::setw(20) << std::left << "P_value";
  caese_file_ << std::setw(25) << std::left << "case_carrier";
  caese_file_ << std::setw(25) << std::left << "case_non_carrier";
  caese_file_ << std::setw(25) << std::left << "control_carrier";
  caese_file_ << std::setw(25) << std::left << "control_non_carrier";
  caese_file_ << std::endl;

  if (tp.gene_list) {
    report_caese(res, tp);
  }
}

auto CAESEReporter::report(std::vector<CAESETask> &resv, TaskParams &tp)
    -> void {
  report_caese(resv, tp);
}

auto CAESEReporter::report_caese(std::vector<CAESETask> &resv, TaskParams &tp)
    -> void {
  for (auto &ct : resv) {
    for (auto &cr : ct.results) {
      std::sort(cr.second.permuted.begin(), cr.second.permuted.end());
      // TODO: interpolate when it isn't an integer
      int lo = cr.second.permuted.size() * 0.025;
      int hi = cr.second.permuted.size() * 0.975;
      caese_file_ << std::setw(25) << std::left << cr.second.gene << " ";
      caese_file_ << std::setw(20) << std::left << cr.second.transcript;
      caese_file_ << std::setw(20) << std::left << cr.second.original;
      caese_file_ << std::setw(20) << std::left << cr.second.permuted[lo];
      caese_file_ << std::setw(20) << std::left << cr.second.permuted[hi];
      caese_file_ << std::setw(20) << std::left << cr.second.or_p;
      caese_file_ << std::setw(25) << std::left << cr.second.case_alt;
      caese_file_ << std::setw(25) << std::left << cr.second.case_ref;
      caese_file_ << std::setw(25) << std::left << cr.second.cont_alt;
      caese_file_ << std::setw(25) << std::left << cr.second.cont_ref;
      caese_file_ << std::endl;
    }
  }
}

auto CAESEReporter::sync_write_caese(std::map<std::string, Result> &crv)
    -> void {
  std::unique_lock<std::mutex> lock(lock_); // Acquire lock

  for (auto &cr : crv) {
    std::sort(cr.second.permuted.begin(), cr.second.permuted.end());
    // TODO: interpolate when it isn't an integer
    int lo = cr.second.permuted.size() * 0.025;
    int hi = cr.second.permuted.size() * 0.975;
    caese_file_ << std::setw(25) << std::left << cr.second.gene << " ";
    caese_file_ << std::setw(20) << std::left << cr.second.transcript;
    caese_file_ << std::setw(20) << std::left << cr.second.original;
    caese_file_ << std::setw(20) << std::left << cr.second.permuted[lo];
    caese_file_ << std::setw(20) << std::left << cr.second.permuted[hi];
    caese_file_ << std::setw(20) << std::left << cr.second.or_p;
    caese_file_ << std::setw(25) << std::left << cr.second.case_alt;
    caese_file_ << std::setw(25) << std::left << cr.second.case_ref;
    caese_file_ << std::setw(25) << std::left << cr.second.cont_alt;
    caese_file_ << std::setw(25) << std::left << cr.second.cont_ref;
    caese_file_ << std::endl;
  }
  lock.unlock();
}

auto CAESEReporter::set_ncases(int ncases) -> void { ncases_ = ncases; }

auto CAESEReporter::set_ncontrols(int ncontrols) -> void {
  ncontrols_ = ncontrols;
}

auto CAESEReporter::cleanup(TaskParams &tp) -> void { return; }
