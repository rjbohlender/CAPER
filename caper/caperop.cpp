//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#include "caperop.hpp"
#include "../utility/math.hpp"

CAPEROp::CAPEROp(CAPERTask &ct, std::shared_ptr<Reporter> reporter, double seed,
                 bool verbose)
    : done_(false), verbose_(verbose), caperTask(ct),
      reporter_(std::move(reporter)) {}

CAPEROp::CAPEROp(CAPERTask &&ct, std::shared_ptr<Reporter> reporter,
                 double seed, bool verbose)
    : done_(false), verbose_(verbose), caperTask(std::move(ct)),
      reporter_(std::move(reporter)) {}

auto CAPEROp::run() -> void {
  Stage stage = caperTask.stage;

  if (stage == Stage::Stage1) {
    op();
  } else if (stage == Stage::Done) {
    caperTask.methods.clear(caperTask.gene.get_transcripts());
    this->done_ = true;
    return;
  } else {
    throw(std::runtime_error("Incorrect stage in CAPEROp::run()"));
  }
}

auto CAPEROp::finish() -> void {
  caperTask.calc_multitranscript_pvalues();
  // Finalize result calculations and free memory
  caperTask.cleanup();

  // Report results for non-gene_list run
  if (!caperTask.tp.gene_list) {
    reporter_->sync_write_simple(caperTask.results, caperTask.tp);
    reporter_->sync_write_detail(caperTask.gene.get_detail(),
                                 caperTask.gene.testable);
    reporter_->sync_write_vaast(caperTask, caperTask.tp);
  }
}

auto CAPEROp::op() -> void {
  auto &gene = caperTask.gene;
  auto &cov = caperTask.get_cov();

  // Set original value
  for (auto &[transcript, res] : caperTask.results) {
    if (!gene.is_polymorphic(transcript)) {
      res.done = true;
      continue;
    }
    if (!res.is_set) {
      res.case_alt = arma::accu(gene.genotypes[transcript].t() *
                                cov.get_original_phenotypes());
      res.case_ref =
          2 * arma::accu(cov.get_original_phenotypes()) - res.case_alt;
      res.cont_alt = arma::accu(gene.genotypes[transcript].t() *
                                (1. - cov.get_original_phenotypes()));
      res.cont_ref =
          2 * arma::accu(1. - cov.get_original_phenotypes()) - res.cont_alt;
      res.original = caperTask.methods.call(
          gene, cov, cov.get_original_phenotypes(), transcript, true);
      res.is_set = true;
    }
  }

  int iter = 0;
  while (iter < caperTask.npermutations) {
    int transcript_no = -1;

    for (auto &[transcript, res] : caperTask.results) {
      transcript_no++;

      // Skip transcripts with no variants
      if (!gene.is_polymorphic(transcript)) {
        continue;
      }

      arma::vec phenotypes = arma::conv_to<arma::vec>::from(
          caperTask.get_permutations()[caperTask.offset + iter]);
#ifndef NDEBUG
      if (arma::accu(phenotypes) != arma::accu(cov.get_original_phenotypes())) {
        std::cerr
            << "vec count: "
            << std::accumulate(
                   caperTask.get_permutations()[caperTask.offset + iter]
                       .begin(),
                   caperTask.get_permutations()[caperTask.offset + iter].end(),
                   0)
            << std::endl;
        std::cerr << "Permuted count: " << arma::accu(phenotypes) << std::endl;
        std::cerr << "Original count: "
                  << arma::accu(caperTask.cov->get_original_phenotypes())
                  << std::endl;
        throw(std::runtime_error("Failed to properly generate permutation."));
      }
#endif

      double perm_val =
          caperTask.methods.call(gene, cov, phenotypes, transcript, false);

      res.permuted.push_back(perm_val);

      // Update total number of permutations
      res.permutations++;

      check_perm(caperTask, caperTask.tp, perm_val, caperTask.success_threshold,
                 res, caperTask.termination);
    }
    // Stop iterating if everything is done
    if (std::all_of(caperTask.results.cbegin(), caperTask.results.cend(),
                    [](const auto &v) { return v.second.done; })) {
      break;
    }
    iter++;
  }

  int ts_no = 0;
  for (auto &[transcript, res] : caperTask.results) {
    // For runs with 0 nperm
    check_done(caperTask.tp, caperTask.success_threshold, res,
               caperTask.termination);
    double empirical;
    double midp;
    // If we stopped early, use the geometric correction, otherwise calculate
    // for binomial p-hat
    if (caperTask.tp.max_perms) { // We're looping multiple times
      if (res.permutations < *caperTask.tp.max_perms) {
        if (res.permutations <= 1) {
          empirical = binomial_estimate(res.successes, res.permutations);
        } else {
          empirical = geometric_p(res.successes, res.permutations);
        }
      } else {
        empirical = binomial_estimate(res.successes, res.permutations);
      }
      midp = (0.5 + res.mid_successes) / (1. + res.permutations);
    } else { // We're looping a single time
      if (res.permutations < caperTask.tp.nperm) {
        if (res.permutations <= 1) {
          empirical = binomial_estimate(res.successes, res.permutations);
        } else {
          empirical = geometric_p(res.successes, res.permutations);
        }
      } else {
        empirical = binomial_estimate(res.successes, res.permutations);
      }
      midp = (0.5 + res.mid_successes) / (1. + res.permutations);
    }

    // Success on every iteration
    if (empirical > 1) {
      empirical = 1;
    }

    if (midp > 1) {
      midp = 1;
    }

    if (ts_no == 0) {
      done_ = res.done;
    } else {
      done_ &= res.done;
    }
    res.empirical_p = empirical;
    res.empirical_midp = midp;
    res.update_ci();

    double n = 2. * gene.get_samples().size();
    const arma::sp_mat column_sums =
        arma::sum(gene.genotypes[transcript], 1);
    double nmac = static_cast<double>(column_sums.n_nonzero);
    double nmaj = n - nmac;

    res.calc_exact_p(nmac, nmaj);
    ts_no++;
  }
  // True only when all are finished, or we've exhausted permutations
  if (done_) {
    if (verbose_) {
      for (const auto &[ts, result] : caperTask.results) {
        std::cerr << gene.gene_name << "\t" << ts << "\t";
        std::cerr << std::defaultfloat << std::setprecision(6)
                  << result.original << std::endl;
      }
    }
    caperTask.stage = Stage::Done;
  } else {
    caperTask.stage = Stage::Stage1;
  }
}

auto CAPEROp::check_perm(const CAPERTask &ct, const TaskParams &tp,
                         double perm_val, long success_threshold, Result &res,
                         unsigned long termination) -> void {
  // Some methods return a pvalue, so we need to reverse the success inequality
  if (tp.analytic) {
    if (perm_val <= res.original) {
      res.successes++;
      if (tp.pthresh) {
        if (perm_val == res.original) {
          res.mid_successes += 0.5;
        } else {
          res.mid_successes++;
        }

        double lower, upper;
        std::tie(lower, upper) = poisson_ci(res.successes, res.permutations);

        // Ensure a minimum number of permutations
        if (res.permutations > 10) {
          // Stop when the lower bound on the 95% confidence interval is greater
          // than the threshold given
          res.done = lower > *tp.pthresh;
        }
      } else if (res.successes < success_threshold) {
        if (perm_val == res.original) {
          res.mid_successes += 0.5;
        } else {
          res.mid_successes++;
        }
      } else {
        if (perm_val == res.original) {
          res.mid_successes += 0.5;
        } else {
          res.mid_successes++;
        }
        res.done = true;
      }
    }
  } else {
    if (perm_val >= res.original) {
      res.successes++;
      if (tp.pthresh) {
        if (perm_val == res.original) {
          res.mid_successes += 0.5;
        } else {
          res.mid_successes++;
        }

        double lower, upper;
        std::tie(lower, upper) = poisson_ci(res.successes, res.permutations);

        // Ensure a minimum number of permutations
        if (res.permutations > 10) {
          // Stop when the lower bound on the 95% confidence interval is greater
          // than the threshold given
          res.done = lower > *tp.pthresh;
        }
      } else if (res.successes < success_threshold) {
        if (perm_val == res.original) {
          res.mid_successes += 0.5;
        } else {
          res.mid_successes++;
        }
      } else {
        if (perm_val == res.original) {
          res.mid_successes += 0.5;
        } else {
          res.mid_successes++;
        }
        res.done = true;
      }
    }
  }
  if (tp.max_perms) {
    if (tp.gene_list) {
      res.done |= res.permutations >= termination;
      // res.done |= (res.permutations >= (*tp.max_perms / (tp.nthreads - 1)));
    } else {
      res.done |= res.permutations >= *tp.max_perms;
    }
  } else {
    if (tp.gene_list) {
      res.done |= (res.permutations >= ct.npermutations);
    } else {
      res.done |= res.permutations >= tp.nperm;
    }
  }
}

auto CAPEROp::check_done(const TaskParams &tp, long success_threshold,
                         Result &res, unsigned long termination) -> void {
  // Some methods return a pvalue, so we need to reverse the success inequality
  if (tp.analytic) {
    if (tp.pthresh) {
      double lower, upper;
      std::tie(lower, upper) = poisson_ci(res.successes, res.permutations);

      // Ensure a minimum number of permutations
      if (res.permutations > 10) {
        // Stop when the lower bound on the 95% confidence interval is greater
        // than the threshold given
        res.done = lower > *tp.pthresh;
      }
    } else if (res.successes >= success_threshold) {
      res.done = true;
    }
  } else {
    if (tp.pthresh) {
      double lower, upper;
      std::tie(lower, upper) = poisson_ci(res.successes, res.permutations);

      // Ensure a minimum number of permutations
      if (res.permutations > 10) {
        // Stop when the lower bound on the 95% confidence interval is greater
        // than the threshold given
        res.done = lower > *tp.pthresh;
      }
    } else if (res.successes >= success_threshold) {
      res.done = true;
    }
  }
  if (tp.max_perms) {
    if (tp.gene_list) {
      res.done |= res.permutations >= termination;
      const auto worker_threads = tp.nthreads > 1 ? tp.nthreads - 1 : 0;
      if (worker_threads > 0) {
        res.done |=
            (res.permutations >= (*tp.max_perms / worker_threads));
      } else {
        res.done |= res.permutations >= *tp.max_perms;
      }
    } else {
      res.done |= res.permutations >= *tp.max_perms;
    }
  } else {
    if (tp.gene_list) {
      res.done |= res.permutations >= termination;
    } else {
      res.done |= res.permutations >= tp.nperm;
    }
  }
}
