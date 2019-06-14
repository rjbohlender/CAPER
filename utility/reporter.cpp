//
// Created by Bohlender,Ryan James on 10/11/18.
//

#include "reporter.hpp"
#include "filesystem.hpp"

#include <boost/format.hpp>

const std::set<std::string> Reporter::pvalue_methods_ {
  "SKAT",
  "SKATO",
  "CMC",
  "RVT1",
  "RVT2"
};

Reporter::Reporter(TaskParams &tp)
: method_(tp.method), gene_list_(tp.gene_list), testable_(tp.testable), ncases_(0), ncontrols_(0) {
  std::string header;

  if(!check_directory_exists(tp.output_path)) {
    throw(std::runtime_error("Output path is invalid."));
  }

  // If using VAAST, output VAAST file
  if(tp.method == "VAAST") {
    vaast_path_ss << tp.output_path << "/" << tp.method;
	if(tp.group_size > 0)
	  vaast_path_ss << boost::format(".g%1$d") % tp.group_size;
	if(tp.testable)
	  vaast_path_ss << ".testable";
	if(tp.biallelic)
	  vaast_path_ss << ".biallelic";
	if (tp.range_start && tp.range_end) {
	  vaast_path_ss << "." << *tp.range_start;
	  vaast_path_ss << "." << *tp.range_end;
	}
	vaast_path_ss << ".vaast";
	vaast_file_ = std::ofstream(vaast_path_ss.str());
  }

  simple_path_ss << tp.output_path << "/" << tp.method;
  if(tp.group_size > 0)
    simple_path_ss << boost::format(".g%1$d") % tp.group_size;
  if(tp.testable)
    simple_path_ss << ".testable";
  if(tp.biallelic)
    simple_path_ss << ".biallelic";
  if (tp.range_start && tp.range_end) {
	simple_path_ss << "." << *tp.range_start;
	simple_path_ss << "." << *tp.range_end;
  }
  simple_path_ss << ".simple";
  simple_path_tmp_ss << simple_path_ss.str() << ".tmp";

  simple_file_tmp_ = std::ofstream(simple_path_tmp_ss.str());
  if(!simple_file_tmp_.good()) {
    throw(std::runtime_error("Simple file failed to open for writing.\n"));
  }
  simple_file_tmp_ << std::setw(20) << std::left << "Gene";
  simple_file_tmp_ << std::setw(20) << "Transcript";
  simple_file_tmp_ << std::setw(20) << "Test_Statistic";
  simple_file_tmp_ << std::setw(20) << "Empirical_P";
  simple_file_tmp_ << std::setw(20) << "Empirical_MidP";
  simple_file_tmp_ << std::setw(20) << "MGIT_P";
  simple_file_tmp_ << std::setw(20) << "Successes";
  simple_file_tmp_ << std::setw(20) << "MGIT_Successes";
  simple_file_tmp_ << std::setw(20) << "Permutations" << std::endl;

  detail_path_ss << tp.output_path << "/" << tp.method;
  if(tp.group_size > 0)
    detail_path_ss << boost::format(".g%1$d") % tp.group_size;
  if(tp.testable)
    detail_path_ss << ".testable";
  if(tp.biallelic)
	detail_path_ss << ".biallelic";
  if (tp.range_start && tp.range_end) {
	detail_path_ss << "." << *tp.range_start;
	detail_path_ss << "." << *tp.range_end;
  }
  detail_path_ss << ".detail";

  detail_file_ = std::ofstream(detail_path_ss.str());
  if(!detail_file_.good()) {
    throw(std::runtime_error("Detail file failed to open for writing.\n"));
  }

  if(tp.linear) {
    header =
        "Gene\tTranscripts\tVariant\tScore\tOR\tOR_SE\tOR_P\tAF";
  } else {
    header =
        "Gene\tTranscripts\tVariant\tScore\tOR\tOR_SE\tOR_P\tAF\tcase_ref\tcase_alt\tcontrol_ref\tcontrol_alt\tcase_list\tcontrol_list";
  }
  detail_file_ << header << std::endl;
}

Reporter::Reporter(std::vector<CARVATask> &res, TaskParams &tp)
: method_(tp.method), gene_list_(tp.gene_list), testable_(tp.testable), ncases_(0), ncontrols_(0) {
  if(!check_directory_exists(tp.output_path)) {
    throw(std::runtime_error("Output path is invalid."));
  }

  // Normal execution
  simple_path_ss << tp.output_path << "/" << tp.method;
  if (tp.group_size > 0)
	simple_path_ss << boost::format(".g%1$d") % tp.group_size;
  if (tp.testable)
	simple_path_ss << ".testable";
  if (tp.biallelic)
	simple_path_ss << ".biallelic";
  if (tp.range_start && tp.range_end) {
	simple_path_ss << "." << *tp.range_start;
	simple_path_ss << "." << *tp.range_end;
  }
  simple_path_ss << ".simple";

  detail_path_ss << tp.output_path << "/" << tp.method;
  if (tp.group_size > 0)
	detail_path_ss << boost::format(".g%1$d") % tp.group_size;
  if (tp.testable)
	detail_path_ss << ".testable";
  if (tp.biallelic)
	detail_path_ss << ".biallelic";
  if (tp.range_start && tp.range_end) {
	detail_path_ss << "." << *tp.range_start;
	detail_path_ss << "." << *tp.range_end;
  }
  detail_path_ss << ".detail";

  simple_file_tmp_ = std::ofstream(simple_path_ss.str());
  detail_file_ = std::ofstream(detail_path_ss.str());

  if (!simple_file_tmp_.good()) {
	throw (std::runtime_error("Simple file failed to open for writing.\n"));
  }
  if (!detail_file_.good()) {
	throw (std::runtime_error("Detail file failed to open for writing.\n"));
  }
  // Extract results
  extract_results(res, tp);


  // Write output
  if (tp.gene_list) {
	report_simple(tp);
	if (!tp.nodetail) {
	  report_detail(res, tp);
	}
  }
}

auto Reporter::report(std::vector<CARVATask> &res, TaskParams &tp) -> void {
  // Extract results
  extract_results(res, tp);

  // Write output
  if (tp.gene_list) {
	report_simple(tp);
	if (!tp.nodetail) {
	  report_detail(res, tp);
	  if (tp.method == "VAAST") {
		report_vaast(res, tp);
	  }
	}
  }
}

auto Reporter::cleanup(TaskParams &tp) -> void {
  sort_simple(tp);
}

auto Reporter::extract_results(std::vector<CARVATask> &tq_results, TaskParams &tp) -> void {
  if(gene_list_) {
    // Combine results
    std::map<std::string, std::map<std::string, Result>> results;
    for(auto &v : tq_results) {
      for(auto &r : v.results) {
        if (results.find(r.second.gene) != results.end()) {
          if(results[r.second.gene].find(r.second.transcript) != results[r.second.gene].end()) {
            results[r.second.gene][r.second.transcript].combine(r.second);
          } else {
            results[r.second.gene][r.second.transcript] = r.second;
          }
        } else {
          results[r.second.gene][r.second.transcript] = r.second;
		  details_.push_back(v.get_gene().get_detail());
        }
      }
    }
    // Recalculate MGIT
    recalculate_mgit(results);
    // Extract combined results
    for(auto &v : results) {
      for(auto &r : v.second) {
        results_.emplace_back(std::move(r.second));
      }
    }
  } else {
    for(auto &v : tq_results) {
      for(auto &r : v.results) {
        results_.emplace_back(std::move(r.second));
        details_.push_back(v.get_gene().get_detail());
      }
    }
  }

  // Sort extracted results
  if(pvalue_methods_.find(method_) != pvalue_methods_.end() && (method_ != "SKAT" || tp.total_permutations == 0 )) {
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
}

auto Reporter::recalculate_mgit(std::map<std::string, std::map<std::string, Result>> &results) -> void {
  for(auto &g : results) {
    // Skip MGIT if only a single transcript
    if (g.second.size() == 1) {
      for (auto &v : g.second) {
        v.second.mgit_p = v.second.empirical_p;
        v.second.mgit_successes = v.second.successes;
      }
      continue;
    }

    // Shorthand
    unsigned long n = g.second.size();
    arma::uword max_perm = 0;
    arma::uword i, j, k;
    double successes;

    // Get max_perm
    for (const auto &tr : g.second) {
      if (tr.second.permutations > max_perm) {
        max_perm = tr.second.permutations;
      }
    }

    arma::mat mgit_pval_mat = arma::mat(max_perm + 1, n);
    i = 0;
    for(auto &tr : g.second) {
      int m = tr.second.permuted.size();

      tr.second.permuted.push_back(tr.second.original);
      arma::vec permuted = arma::conv_to<arma::vec>::from(tr.second.permuted);
      arma::vec pvals;

      // SKATO and SKAT return pvalues so reverse success criteria
      if (pvalue_methods_.find(method_) != pvalue_methods_.end()) {
        pvals = rank(permuted, "ascend");
      } else {
        pvals = rank(permuted, "descend");
      }

      pvals /= permuted.n_rows;

      try {
        mgit_pval_mat.col(i) = pvals;
      } catch (const std::logic_error &e) {
        std::cerr << "n_row: " << mgit_pval_mat.n_rows << " n_col: " << mgit_pval_mat.n_cols << "\n";
        throw (e);
      }
      i++;
    }
    arma::vec mgit_pval_dist_ = arma::min(mgit_pval_mat, 1);

    double min_p = std::numeric_limits<double>::max();
    for (auto &tr : g.second) {
      if(tr.second.empirical_p < min_p) {
        min_p = tr.second.empirical_p;
      }
    }

    for (auto &tr : g.second) {
      unsigned long m = mgit_pval_dist_.n_rows;  // Total permutations

      successes = arma::find(mgit_pval_dist_ <= min_p).eval().n_rows;

      // Store multi-transcript p-value
      tr.second.mgit_p = (1.0 + successes) / (1.0 + m);
      tr.second.mgit_successes = static_cast<int>(successes);
    }
  }
}

auto Reporter::recalculate_mgit(std::unordered_map<std::string, Result> &results) -> void {
  // Skip MGIT if only a single transcript
  if (results.size() == 1) {
	for (auto &v : results) {
	  v.second.mgit_p = v.second.empirical_p;
	  v.second.mgit_successes = v.second.successes;
	}
	return;
  }

  // Shorthand
  unsigned long n = results.size(); // Number of transcripts
  arma::uword max_perm = 0;
  arma::uword i, j, k;
  double successes;

  // Get max_perm
  for (const auto &tr : results) {
	if (tr.second.permutations > max_perm) {
	  max_perm = tr.second.permutations;
	}
  }

  arma::mat mgit_pval_mat = arma::mat(max_perm + 1, n);
  i = 0;
  for(auto &tr : results) {
	int m = tr.second.permuted.size();

	tr.second.permuted.push_back(tr.second.original);
	arma::vec permuted = arma::conv_to<arma::vec>::from(tr.second.permuted);
	arma::vec pvals;

	// SKATO and SKAT return pvalues so reverse success criteria
	if (pvalue_methods_.find(method_) != pvalue_methods_.end()) {
	  pvals = rank(permuted, "ascend");
	} else {
	  pvals = rank(permuted, "descend");
	}

	pvals /= permuted.n_rows;

	try {
	  mgit_pval_mat.col(i) = pvals;
	} catch (const std::logic_error &e) {
	  std::cerr << "n_row: " << mgit_pval_mat.n_rows << " n_col: " << mgit_pval_mat.n_cols << "\n";
	  throw (e);
	}
	i++;
  }
  arma::vec mgit_pval_dist_ = arma::min(mgit_pval_mat, 1);

  double min_p = std::numeric_limits<double>::max();
  for (auto &tr : results) {
	if(tr.second.empirical_p < min_p) {
	  min_p = tr.second.empirical_p;
	}
  }

  for (auto &tr : results) {
	unsigned long m = mgit_pval_dist_.n_rows;  // Total permutations

	successes = arma::find(mgit_pval_dist_ <= min_p).eval().n_rows;

	// Store multi-transcript p-value
	tr.second.mgit_p = (1.0 + successes) / (1.0 + m);
	tr.second.mgit_successes = static_cast<int>(successes);
  }
}

auto Reporter::report_simple(TaskParams &tp) -> void {
  // Print header and formatted results
  double permutation_mean = 0;
  double permutation_variance = 0;

  for (auto &v : results_) {
	if (v.successes < tp.success_threshold) {
	  if(std::find(unfinished_.begin(), unfinished_.end(), v.gene) == unfinished_.end())
		unfinished_.push_back(v.gene);
	}
	permutation_mean += v.permutations;
  }
  permutation_mean /= results_.size();
  for (auto &v : results_) {
	permutation_variance += std::pow(v.permutations - permutation_mean, 2) / results_.size();
  }
  // Report results
  int rank = 1;
  for(auto &v : results_) {
    v.set_rank(rank);
    write_to_stream(simple_file_tmp_, v);
    rank++;
  }

  simple_file_tmp_.flush();
  simple_file_tmp_.close();

  std::cerr << "Permutation mean: " << permutation_mean << std::endl;
  std::cerr << "Permutation sd: " << std::sqrt(permutation_variance) << std::endl;
  std::cerr << "Transcripts submitted: " << results_.size() << std::endl;
}

auto Reporter::report_detail(std::vector<CARVATask> &res, TaskParams &tp) -> void {
  std::string header;

  detail_file_ << header << std::endl;

  int i = 0; // For each gene
  for (auto &v : details_) {
	detail_file_ << v;
	// Print sample / index map at the end
	if (i == res.size() - 1) {
	  detail_file_ << "## Sample Index Map" << std::endl;
	  int j = 0;
	  for (auto &s : res[0].get_gene().get_samples()) {
		detail_file_ << "#\t" << s << "\t" << j << std::endl;
		j++;
	  }
	}
	i++;
  }
}

auto Reporter::report_vaast(std::vector<CARVATask> &res, TaskParams &tp) -> void {
  // Header information
  vaast_file_ << "## PA_VERSION\t0.0" << std::endl;
  vaast_file_ << "## COMMAND\t" << tp.full_command << std::endl;

  for (auto &ta : res) {
	for(auto &ts : ta.get_gene().get_vaast()) {
	  double z2 = std::pow(1.96, 2);
	  double n = ta.results[ts.first].permutations;
	  double ns = ta.results[ts.first].successes;
	  double nf = ta.results[ts.first].permutations - ta.results[ts.first].successes;
	  double p = (ns + 1.) / (n + 1.);
	  double sp = (p + z2 / (2 * (n + 1.))) / (1. + z2 / (n + 1.));
	  double ci = 1.96 / (1. + z2 / n) * std::sqrt((p * (1 - p) / (n + 1.) + z2 / (4 * n * n)));
	  vaast_file_ << ts.second;
	  vaast_file_ << "SCORE: " << boost::format("%1$.2f") % ta.results[ts.first].original << std::endl;
	  vaast_file_ << "genome_permutation_p: " << p << std::endl;
	  vaast_file_ << "genome_permutation_p_ci: " << ((sp - ci > 0) ? sp - ci : 0) << "," << ((sp + ci > 1 ) ? 1 : sp + ci) << std::endl;
	  vaast_file_ << "num_permutations: " << ta.results[ts.first].permutations << std::endl;
	  vaast_file_ << "total_success: " << ta.results[ts.first].successes << std::endl;
	}
  }
}

auto Reporter::write_to_stream(std::ostream &os, Result &res) -> void {
  if(testable_) {
    if(res.testable) {
      os << res;
    }
  } else {
    os << res;
  }
}

auto Reporter::sync_write_simple(std::unordered_map<std::string, Result> &results, TaskParams &tp, bool top_only) -> void {
  std::unique_lock<std::mutex> lock(lock_); // Acquire lock

  recalculate_mgit(results);

  for(const auto &v : results) {
	if (v.second.successes < tp.success_threshold) {
	  if(std::find(unfinished_.begin(), unfinished_.end(), v.second.gene) == unfinished_.end())
		unfinished_.push_back(v.second.gene);
	}
  }
  if(top_only) {
	// Find the most significant result.
	std::unique_ptr<Result> topres;
	for(auto &r : results) {
	  if(topres == nullptr) {
		topres = std::make_unique<Result>(r.second);
	  } else {
		if(topres->empirical_p > r.second.empirical_p) {
		  topres = std::make_unique<Result>(r.second);
		}
	  }
	}
	simple_file_tmp_ << *topres;
    lock_.unlock();
    return;
  }

  if(testable_) {
    for(auto &tr : results) {
      if(tr.second.testable) {
        simple_file_tmp_ << tr.second;
      }
    }
  } else {
	for(auto &tr : results) {
	  simple_file_tmp_ << tr.second;
	}
  }

  lock.unlock();
}

auto Reporter::sync_write_detail(const std::string &d, bool testable) -> void {
  std::unique_lock<std::mutex> lock(lock_); // Acquire lock

  if(testable_) {
    if(testable)
	  detail_file_ << d;
  } else {
    detail_file_ << d;
  }

  detail_file_.flush();
  lock.unlock();
}

auto Reporter::sort_simple(TaskParams &tp) -> void {
  simple_file_ = std::ofstream(simple_path_ss.str());
  if(!simple_file_.good()) {
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
  simple_file_ << std::setw(20) << "Gene";
  simple_file_ << std::setw(20) << "Transcript";
  simple_file_ << std::setw(20) << "Test_Statistic";
  simple_file_ << std::setw(20) << "Empirical_P";
  simple_file_ << std::setw(20) << "Empirical_MidP";
  simple_file_ << std::setw(20) << "MGIT_P";
  simple_file_ << std::setw(20) << "Successes";
  simple_file_ << std::setw(20) << "MGIT_Successes";
  simple_file_ << std::setw(20) << "Permutations" << std::endl;

  std::vector<ResultLine> results;
  unsigned long lineno = 0;

  while(std::getline(ifs, line)) {
    if(lineno == 0) {
      lineno++;
      continue;
    }
    RJBUtil::Splitter<std::string> splitter(line, " \t");

    ResultLine rs = ResultLine {
      .gene = splitter[0],
      .transcript = splitter[1],
      .original = std::stod(splitter[2]),
      .empirical_p = std::stod(splitter[3]),
      .empirical_midp = std::stod(splitter[4]),
      .mgit = std::stod(splitter[5]),
      .successes = std::stoul(splitter[6]),
      .mgit_successes = std::stoul(splitter[7]),
      .permutations = std::stoul(splitter[8])
    };

    results.push_back(rs);

    lineno++;
  }

  ifs.close();

  // Sort results and write them out
  if(tp.method == "SKATO") {
	std::sort(results.begin(), results.end(), [](ResultLine &a, ResultLine &b) { return a.empirical_p < b.empirical_p; });
  } else {
	std::sort(results.begin(), results.end(), [](ResultLine &a, ResultLine &b) { return a.original < b.original; });
  }

  for(const auto &rs : results) {
    simple_file_ << std::setw(20) << std::left <<  rank;
    simple_file_ << std::setw(20) << rs.gene;
    simple_file_ << std::setw(20) << rs.transcript;
    simple_file_ << std::setw(20) << rs.original;
    simple_file_ << std::setw(20) << rs.empirical_p;
    simple_file_ << std::setw(20) << rs.empirical_midp;
	simple_file_ << std::setw(20) << rs.mgit;
    simple_file_ << std::setw(20) << rs.successes;
    simple_file_ << std::setw(20) << rs.mgit_successes;
	simple_file_ << std::setw(20) << rs.permutations << std::endl;

    rank++;
  }

  simple_file_.flush();
  simple_file_.close();

  // Delete tmp file
  std::remove(simple_path_tmp_ss.str().c_str());

  if (!unfinished_.empty() && tp.total_permutations > 0) {
	// Print command to run unfinished
	std::stringstream uf_ss;
	uf_ss << tp.program_path << " ";
	uf_ss << "-i " << tp.genotypes_path << " ";
	uf_ss << "-c " << tp.covariates_path << " ";
	uf_ss << "-o " << tp.output_path << " ";
	uf_ss << "-p " << tp.ped_path << " ";
	if (tp.bed) {
	  uf_ss << "-b " << *tp.bed << " ";
	}
	if (tp.weight) {
	  uf_ss << "-w " << *tp.weight << " ";
	}
	uf_ss << "-1 0 "; // Skip stage 1
	uf_ss << "-2 " << tp.total_permutations * 10 << " ";
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
	if (tp.mac != 250) {
	  uf_ss << "--mac " << tp.mac << " ";
	}
	if (tp.nodetail) {
	  uf_ss << "--nodetail ";
	}
	if(tp.score_only_minor) {
	  uf_ss << "--score_only_minor ";
	}
	if(tp.score_only_alternative) {
	  uf_ss << "--score_only_alternative ";
	}

	uf_ss << "-l ";
	for (int i = 0; i < unfinished_.size(); i++) {
	  if (i == unfinished_.size() - 1) {
		uf_ss << unfinished_[i];
	  } else {
		uf_ss << unfinished_[i] << ",";
	  }
	}
	std::cerr << "Some genes did not reach the success threshold. Run the following command to check those genes."
			  << std::endl;
	std::cerr << uf_ss.str() << std::endl;
  }
}

auto Reporter::set_ncases(int ncases) -> void {
  ncases_ = ncases;
}

auto Reporter::set_ncontrols(int ncontrols) -> void {
  ncontrols_ = ncontrols;
}

PowerReporter::PowerReporter(TaskParams &tp) {
  std::stringstream power_path_ss;

  if(tp.method == "VAAST" && tp.group_size > 0) {
	power_path_ss << tp.output_path << "/" << tp.method << ".g" << tp.group_size << ".power";
  } else {
	power_path_ss << tp.output_path << "/" << tp.method << ".power";
  }
  power_file_ = std::ofstream(power_path_ss.str());

// Initial header write
  power_file_ << std::setw(15) << "Gene";
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

  if(tp.method == "VAAST" && tp.group_size > 0) {
	power_path_ss << tp.output_path << "/" << tp.method << ".g" << tp.group_size << ".power";
  } else {
	power_path_ss << tp.output_path << "/" << tp.method << ".power";
  }
  power_file_ = std::ofstream(power_path_ss.str());

// Initial header write
  power_file_ << std::setw(15) << "Gene";
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

auto PowerReporter::report(std::vector<PowerTask> &resv, TaskParams &tp) -> void {
  report_power(resv, tp);
}


auto PowerReporter::report_power(std::vector<PowerTask> &resv, TaskParams &tp) -> void {
  for (const auto &pt : resv) {
    for (const auto &pr : pt.power_res_) {
	  power_file_ << std::setw(20) << pr.gene;
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

  for(const auto &pr : prv) {
	power_file_ << std::setw(20) << pr.gene;
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

auto PowerReporter::set_ncases(int ncases) -> void {
  ncases_ = ncases;
}

auto PowerReporter::set_ncontrols(int ncontrols) -> void {
  ncontrols_ = ncontrols;
}

auto PowerReporter::cleanup(TaskParams &tp) -> void {
  return;
}

CAESEReporter::CAESEReporter(TaskParams &tp) {
  std::stringstream caese_path_ss;

  if(tp.method == "VAAST" && tp.group_size > 0) {
	caese_path_ss << tp.output_path << "/" << tp.method << ".g" << tp.group_size << ".caese";
  } else {
	caese_path_ss << tp.output_path << "/" << tp.method << ".caese";
  }
  caese_file_ = std::ofstream(caese_path_ss.str());

// Initial header write
  caese_file_ << std::setw(20) << std::left << "Gene";
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

  if(tp.method == "VAAST" && tp.group_size > 0) {
	caese_path_ss << tp.output_path << "/" << tp.method << ".g" << tp.group_size << ".caese";
  } else {
	caese_path_ss << tp.output_path << "/" << tp.method << ".caese";
  }
  caese_file_ = std::ofstream(caese_path_ss.str());

// Initial header write
  caese_file_ << std::setw(20) << std::left << "Gene";
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

auto CAESEReporter::report(std::vector<CAESETask> &resv, TaskParams &tp) -> void {
  report_caese(resv, tp);
}

auto CAESEReporter::report_caese(std::vector<CAESETask> &resv, TaskParams &tp) -> void {
  for (auto &ct : resv) {
	for (auto &cr : ct.results) {
	  std::sort(cr.second.permuted.begin(), cr.second.permuted.end());
	  // TODO: interpolate when it isn't an integer
	  int lo = cr.second.permuted.size() * 0.025;
	  int hi = cr.second.permuted.size() * 0.975;
	  caese_file_ << std::setw(20) << std::left << cr.second.gene;
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

auto CAESEReporter::sync_write_caese(std::map<std::string, Result> &crv) -> void {
  std::unique_lock<std::mutex> lock(lock_); // Acquire lock

  for (auto &cr : crv) {
	std::sort(cr.second.permuted.begin(), cr.second.permuted.end());
	// TODO: interpolate when it isn't an integer
	int lo = cr.second.permuted.size() * 0.025;
	int hi = cr.second.permuted.size() * 0.975;
	caese_file_ << std::setw(20) << std::left << cr.second.gene;
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

auto CAESEReporter::set_ncases(int ncases) -> void {
  ncases_ = ncases;
}

auto CAESEReporter::set_ncontrols(int ncontrols) -> void {
  ncontrols_ = ncontrols;
}

auto CAESEReporter::cleanup(TaskParams &tp) -> void {
  return;
}
