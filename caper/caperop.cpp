//
// Created by Bohlender,Ryan James on 2019-05-31.
//

#include "caperop.hpp"
#include "../utility/math.hpp"

CAPEROp::CAPEROp(CAPERTask &ct, std::shared_ptr<Reporter> reporter, double seed,
                 bool verbose)
    : done_(false), verbose_(verbose), carvaTask(ct),
      reporter_(std::move(reporter)) {}

CAPEROp::CAPEROp(CAPERTask &&ct, std::shared_ptr<Reporter> reporter,
                 double seed, bool verbose)
    : done_(false), verbose_(verbose), carvaTask(std::move(ct)),
      reporter_(std::move(reporter)) {}

auto CAPEROp::run() -> void {
  Stage stage = carvaTask.stage;

  if (stage == Stage::Stage1) {
    op();
  } else if (stage == Stage::Done) {
    carvaTask.methods.clear(carvaTask.gene.get_transcripts());
    this->done_ = true;
    return;
  } else {
    throw(std::runtime_error("Incorrect stage in CAPEROp::run()"));
  }
}

auto CAPEROp::finish() -> void {
  carvaTask.calc_multitranscript_pvalues();
  // Finalize result calculations and free memory
  carvaTask.cleanup();

  // Report results for non-gene_list run
  if (!carvaTask.tp.gene_list) {
    reporter_->sync_write_simple(carvaTask.results, carvaTask.tp);
    reporter_->sync_write_detail(carvaTask.gene.get_detail(),
                                 carvaTask.gene.testable);
    reporter_->sync_write_vaast(carvaTask, carvaTask.tp);
  }
}

auto CAPEROp::op() -> void {
  auto &gene = carvaTask.gene;
  auto &cov = carvaTask.get_cov();

  // Set original value
  for (auto &[transcript, res] : carvaTask.results) {
    if (!gene.is_polymorphic(transcript)) {
      res.done = true;
      continue;
    }
    if (!res.is_set) {
      res.case_alt = arma::accu(gene.genotypes[transcript].t() *
                                cov.get_original_phenotypes());
      res.case_ref = 2 * arma::accu(cov.get_original_phenotypes()) - res.case_alt;
      res.cont_alt = arma::accu(gene.genotypes[transcript].t() *
                                (1. - cov.get_original_phenotypes()));
      res.cont_ref =
          2 * arma::accu(1. - cov.get_original_phenotypes()) - res.case_alt;
      res.original = carvaTask.methods.call(
          gene, cov, cov.get_original_phenotypes(), transcript, true);
      res.is_set = true;
    }
  }

  int iter = 0;
  while (iter < carvaTask.npermutations) {
    int transcript_no = -1;

    for (auto &[transcript, res] : carvaTask.results) {
      transcript_no++;

      // Skip transcripts with no variants
      if (!gene.is_polymorphic(transcript)) {
        continue;
      }

      arma::vec phenotypes = arma::conv_to<arma::vec>::from(
          carvaTask.get_permutations()[carvaTask.offset + iter]);
#ifndef NDEBUG
      if (arma::accu(phenotypes) != arma::accu(cov.get_original_phenotypes())) {
        std::cerr
            << "vec count: "
            << std::accumulate(
                   carvaTask.get_permutations()[carvaTask.offset + iter]
                       .begin(),
                   carvaTask.get_permutations()[carvaTask.offset + iter].end(),
                   0)
            << std::endl;
        std::cerr << "Permuted count: " << arma::accu(phenotypes) << std::endl;
        std::cerr << "Original count: "
                  << arma::accu(carvaTask.cov->get_original_phenotypes())
                  << std::endl;
        throw(std::runtime_error("Failed to properly generate permutation."));
      }
#endif

      double perm_val =
          carvaTask.methods.call(gene, cov, phenotypes, transcript, false);

      res.permuted.push_back(perm_val);

      // Update total number of permutations
      res.permutations++;

      check_perm(carvaTask.tp, perm_val, carvaTask.success_threshold, transcript, res);
    }
    // Stop iterating if everything is done
    if (std::all_of(carvaTask.results.cbegin(), carvaTask.results.cend(),
                    [](const auto &v) { return v.second.done; })) {
      break;
    }
    iter++;
  }

  int ts_no = 0;
  for (auto &[transcript, res] : carvaTask.results) {
    double empirical;
    double midp;
    // If we stopped early, use the geometric correction, otherwise calculate
    // for binomial p-hat
    if (carvaTask.tp.max_perms) { // We're looping multiple times
      if (res.permutations < *carvaTask.tp.max_perms) {
        empirical = geometric_p(res.successes, res.permutations);
        midp = geometric_p(res.mid_successes,
                           static_cast<double>(res.permutations));
      } else {
        empirical = (1. + res.successes) / (1. + res.permutations);
        midp = (1. + res.mid_successes) / (1. + res.permutations);
      }
    } else { // We're looping a single time
      if (res.permutations < carvaTask.tp.nperm) {
        empirical = geometric_p(res.successes, res.permutations);
        midp = geometric_p(res.mid_successes,
                           static_cast<double>(res.permutations));
      } else {
        empirical = (1. + res.successes) / (1. + res.permutations);
        midp = (1. + res.mid_successes) / (1. + res.permutations);
      }
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
    if (carvaTask.tp.max_perms) {
      if (carvaTask.tp.gene_list) {
        if (res.permutations >= carvaTask.termination) {
          done_ = true;
        }
      } else {
        if (res.permutations >= *carvaTask.tp.max_perms) {
          done_ = true;
        }
      }
    } else {
      if (carvaTask.tp.gene_list) {
        if (res.permutations >= carvaTask.termination) {
          done_ = true;
        }
      } else {
        if (res.permutations >= carvaTask.tp.nperm) {
          done_ = true;
        }
      }
    }
    res.empirical_p = empirical;
    res.empirical_midp = midp;
    res.update_ci();

    double n = 2. * gene.get_samples().size();
    double nmac =
        arma::accu(arma::sum(arma::mat(gene.genotypes[transcript]), 1) > 0);
    double nmaj = n - nmac;

    res.calc_exact_p(nmac, nmaj);
    ts_no++;
  }
  if (done_) { // True only when all are finished.
    if (verbose_) {
      for (const auto &v : carvaTask.results) {
        std::cerr << gene.gene_name << "\t" << v.second.transcript << "\t";
        std::cerr << std::defaultfloat << std::setprecision(6)
                  << v.second.original << std::endl;
      }
    }
    carvaTask.stage = Stage::Done;
  } else {
    carvaTask.stage = Stage::Stage1;
  }
}

auto CAPEROp::check_perm(const TaskParams &tp, double perm_val,
                         long success_threshold,
                         const std::string &ts, Result &res) -> void {
  // Some methods return a pvalue so we need to reverse the success inequality
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
}
