//
// Created by Bohlender,Ryan James on 10/11/18.
//

#include "reporter.hpp"

#include <boost/format.hpp>

const std::set<std::string> Reporter::pvalue_methods_ {
  "SKAT",
  "SKATO",
  "CMC"
};

Reporter::Reporter(TaskQueue &tq, TaskParams &tp)
: method_(tp.method), gene_list_(tp.gene_list), testable_(tp.testable) {
  // Extract results
  extract_results(tq.get_results());

  // Write output
  report_simple(tp);
  if(!tp.nodetail)
	report_detail(tq, tp);
}

auto Reporter::extract_results(std::vector<TaskArgs> &tq_results) -> void {
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
  if(pvalue_methods_.find(method_) != pvalue_methods_.end()) {
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
        return a.original < b.original;
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

    for (auto &tr : g.second) {
      unsigned long m = mgit_pval_dist_.n_rows;  // Total permutations

      successes = arma::find(mgit_pval_dist_ <= tr.second.empirical_p).eval().n_rows;

      // Store multi-transcript p-value
      tr.second.mgit_p = (1.0 + successes) / (1.0 + m);
      tr.second.mgit_successes = static_cast<int>(successes);
    }
  }
}

auto Reporter::report_simple(TaskParams &tp) -> void {
  std::stringstream simple_path_ss;
  simple_path_ss << tp.output_path << "/" << tp.method;
  if(tp.group_size > 0)
    simple_path_ss << boost::format(".g%1$d") % tp.group_size;
  if(tp.testable)
    simple_path_ss << ".testable";
  simple_path_ss << ".simple";
  std::ofstream simple_ofs(simple_path_ss.str());

  // Holds unfinished genes
  std::vector<std::string> unfinished;

  simple_ofs << std::setw(20) << std::left << "Rank";
  simple_ofs << std::setw(20) << "Gene";
  simple_ofs << std::setw(20) << "Transcript";
  simple_ofs << std::setw(20) << "Original";
  simple_ofs << std::setw(20) << "Empirical_P";
  simple_ofs << std::setw(20) << "Empirical_MidP";
  simple_ofs << std::setw(20) << "Successes";
  simple_ofs << std::setw(20) << "Permutations";
  simple_ofs << std::setw(20) << "MGIT";
  simple_ofs << std::setw(20) << "MGIT_Successes";
  simple_ofs << std::setw(20) << "OddsRatio" << std::endl;
  // Print header and formatted results
  double permutation_mean = 0;
  double permutation_variance = 0;

  bool unfin = false;
  for (auto &v : results_) {
	if (v.successes < tp.success_threshold && !unfin) {
	  unfinished.push_back(results_.back().gene);
	  unfin = true;
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
    write_to_stream(simple_ofs, v);
    rank++;
  }

  simple_ofs.flush();
  simple_ofs.close();

  std::cerr << "Permutation mean: " << permutation_mean << std::endl;
  std::cerr << "Permutation sd: " << std::sqrt(permutation_variance) << std::endl;
  std::cerr << "Transcripts submitted: " << results_.size() << std::endl;

  if (!unfinished.empty()) {
    // Print command to run unfinished
    std::stringstream uf_ss;
    uf_ss << tp.program_path << " ";
    uf_ss << "-g " << tp.genotypes_path << " ";
    uf_ss << "-c " << tp.covariates_path << " ";
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
    if (tp.adjust) {
      uf_ss << "-n ";
    }
    if (!tp.verbose) {
      uf_ss << "-q ";
    }
    if (tp.a != 1 || tp.b != 25) {
      uf_ss << "--beta_weights " << tp.a << "," << tp.b << " ";
    }
    if (tp.group_size > 0) {
      uf_ss << " -g " << tp.group_size;
    }

    uf_ss << "-l ";
    for (int i = 0; i < unfinished.size(); i++) {
      if (i == unfinished.size() - 1) {
        uf_ss << unfinished[i];
      } else {
        uf_ss << unfinished[i] << ",";
      }
    }
    std::cerr << "Some genes did not reach the success threshold. Run the following command to check those genes."
              << std::endl;
    std::cerr << uf_ss.str() << std::endl;
  }
}

auto Reporter::report_detail(TaskQueue &tq, TaskParams &tp) -> void {
  std::stringstream detail_path_ss;
  detail_path_ss << tp.output_path << "/" << tp.method;
  if(tp.group_size > 0)
    detail_path_ss << boost::format(".g%1$d") % tp.group_size;
  if(tp.testable)
    detail_path_ss << ".testable";
  detail_path_ss << ".detail";

  std::ofstream detail(detail_path_ss.str());

  std::string header =
      "#Gene\tTranscripts\tVariant\tScore\tOR\tOR_SE\tOR_P\tAF\tcase_ref\tcase_alt\tcontrol_ref\tcontrol_alt\tcase_list\tcontrol_list";
  detail << header << std::endl;

  int i = 0; // For each gene
  for (auto &v : details_) {
	detail << v;
	// Print sample / index map at the end
	if (i == tq.get_results().size() - 1) {
	  detail << "## Sample Index Map" << std::endl;
	  int j = 0;
	  for (auto &s : tq.get_results()[0].get_gene().get_samples()) {
		detail << "#\t" << s << "\t" << j << std::endl;
		j++;
	  }
	}
	i++;
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

