//
// Created by Bohlender,Ryan James on 8/22/18.
//

#include "taskqueue.hpp"

TaskQueue::TaskQueue(size_t thread_cnt)
	: threads_(thread_cnt),
	  ntasks_(0),
	  quit_(false),
	  gen_(rd_()),
	  verbose_(false),
	  nthreads_(thread_cnt) {
  for (size_t i = 0; i < threads_.size(); i++) {
	threads_[i] = std::thread(
		std::bind(&TaskQueue::thread_handler, this));
  }
}

TaskQueue::TaskQueue(size_t thread_cnt, bool verbose)
	: threads_(thread_cnt),
	  ntasks_(0),
	  quit_(false),
	  gen_(rd_()),
	  verbose_(verbose),
	  nthreads_(thread_cnt) {
  for (size_t i = 0; i < threads_.size(); i++) {
	threads_[i] = std::thread(
		std::bind(&TaskQueue::thread_handler, this));
  }
}

TaskQueue::~TaskQueue() {
  quit_ = true;
  cv_.notify_all();

  // Wait for threads to finish
  for (size_t i = 0; i < threads_.size(); i++) {
	if (threads_[i].joinable()) {
	  threads_[i].join();
	}
  }
}

void TaskQueue::join() {
  using namespace std::chrono_literals;
  // Complete jobs
  while (!q_.empty() || ntasks_ > 0) {
	std::this_thread::sleep_for(0.1s);
  }

  quit_ = true;
  cv_.notify_all();

  for (size_t i = 0; i < threads_.size(); i++) {
	if (threads_[i].joinable()) {
	  threads_[i].join();
	}
  }
}

void TaskQueue::dispatch(TaskArgs &ta) {
  std::unique_lock<std::mutex> lock(lock_);

  q_.push(ta);
  ntasks_++;

  lock.unlock();
  cv_.notify_all();
}

void TaskQueue::dispatch(TaskArgs &&ta) {
  std::unique_lock<std::mutex> lock(lock_);

  q_.emplace(ta);
  ntasks_++;

  lock.unlock();
  cv_.notify_all();
}

bool TaskQueue::empty() {
  return q_.empty();
}

std::vector<TaskArgs> &TaskQueue::get_results() {
  return results_;
}

void TaskQueue::thread_handler() {
  std::unique_lock<std::mutex> lock(lock_);
  do {
	// Wait for data
	cv_.wait(lock, [this] {
	  return q_.size() || quit_;
	});

	// After waiting, we have the lock
	if (q_.size() && !quit_) {
	  auto op = q_.front();
	  q_.pop();

	  // unlock, now that we're done with queue
	  lock.unlock();

	  Stage stage = op.get_stage();

	  if (stage == Stage::Stage1) {
		stage_1(op);
		ntasks_--;
	  } else if (stage == Stage::Stage2) {
		stage_2(op);
		ntasks_--;
	  } else if (stage == Stage::Done) {
		// MGIT
		op.calc_multitranscript_pvalues();
		// Free memory
		op.cleanup();

		// Lock and do results queue
		lock.lock();
		results_.emplace_back(op);
		lock.unlock();
		ntasks_--;
	  }

	  lock.lock();
	}
  } while (!quit_ || ntasks_ > 0);
}

void TaskQueue::stage_1(TaskArgs &ta) {
  // Set original value
  for (auto &v : ta.results) {
	v.second.original = call_method(ta.get_methods(), ta.get_gene(), ta.get_cov(), v.second.transcript, ta.get_tp(), false, true);
  }

  if (verbose_) {
	for (const auto &v : ta.results) {
	  std::cerr << "Stage 1: " << ta.get_gene().get_gene() << "\t" << v.second.transcript << "\t";
	  std::cerr << v.second.original << std::endl;
	}
  }

  int iter = 0;
  while (iter < ta.get_npermutations()) {
	if (!(ta.get_methods().str() == "SKAT") && !(ta.get_methods().str() == "SKATO") && !(ta.get_methods().str() == "BURDEN")) {
	  ta.get_cov().set_phenotype_vector(ta.get_permutations()[iter]);
	}
	int transcript_no = -1;

	for (auto &v : ta.results) {
	  transcript_no++;
	  const std::string &k = v.second.transcript;
	  double perm_val;

	  if (!(ta.get_methods().str() == "SKAT") && !(ta.get_methods().str() == "SKATO") && !(ta.get_methods().str() == "BURDEN")) {
		perm_val = call_method(ta.get_methods(), ta.get_gene(), ta.get_cov(), v.second.transcript, ta.get_tp(), false, false);
	  } else {
		if (transcript_no == 0) {
		  perm_val = call_method(ta.get_methods(), ta.get_gene(), ta.get_cov(), v.second.transcript, ta.get_tp(), true, false);
		} else {
		  perm_val = call_method(ta.get_methods(), ta.get_gene(), ta.get_cov(), v.second.transcript, ta.get_tp(), false, false);
		}
	  }

	  // ta.increment_permuted(v.second.transcript, perm_val);
	  v.second.permuted.push_back(perm_val);

	  check_perm(ta.get_methods().str(), perm_val, ta.success_threshold, v);

	  // Update total number of permutations
	  v.second.permutations++;

	  // Track when we reached success threshold
	  if (v.second.successes == ta.success_threshold && v.second.min_success_at < 0) {
		v.second.min_success_at = v.second.permutations;
	  }
	}
	// Stop iterating if everything is done
	if (std::all_of(ta.results.cbegin(), ta.results.cend(), [&](const auto &v) { return v.second.done; }))
	  break;
	iter++;
  }

  if (std::any_of(ta.results.cbegin(), ta.results.cend(), [&](const auto &v) { return !v.second.done; })
	  && ta.get_remaining() > 0) {
	ta.set_stage(Stage::Stage2);
  } else {
	for (auto &v : ta.results) {
	  double empirical;
	  double midp;
	  if (v.second.min_success_at > 0 && v.second.successes == ta.success_threshold + 1) {
		std::uniform_int_distribution<> dis(v.second.min_success_at, v.second.permutations);

		empirical = v.second.successes / (1. + dis(gen_));
		midp = v.second.mid_successes / (1. + dis(gen_));
	  } else {
		empirical = (1. + v.second.successes) / (1. + v.second.permutations);
		midp = (1. + v.second.mid_successes) / (1. + v.second.permutations);
	  }

	  // Success on every iteration
	  if (empirical > 1) {
		empirical = 1;
	  }

	  if (midp > 1) {
		midp = 1;
	  }

	  ta.set_stage(Stage::Done);
	  v.second.empirical_p = empirical;
	  v.second.empirical_midp = midp;
	}
  }
  dispatch(ta);
}

void TaskQueue::stage_2(TaskArgs &ta) {
  // Declare maps
  std::map<std::string, std::vector<int32_t>> mac_case_count;
  std::map<std::string, std::vector<int32_t>> maj_case_count;
  std::map<std::string, arma::colvec> mac_odds;
  std::map<std::string, arma::colvec> maj_odds;
  std::map<std::string, arma::uvec> mac_indices;
  std::map<std::string, arma::uvec> maj_indices;

  double perm_val;
  int transcript_no;
  // For permutation set output
  std::ofstream pset_ofs;


  // Setup
  for (auto &v : ta.results) {
	const std::string &k = v.second.transcript;

	if (std::isnan(v.second.original)) {
	  v.second.original = call_method(ta.get_methods(), ta.get_gene(), ta.get_cov(), k, ta.get_tp(), false, true);
	}
	// Minor allele carrier indices
	mac_indices[k] = arma::find(arma::sum(ta.get_gene().get_matrix(k), 1) > 0);
	maj_indices[k] = arma::find(arma::sum(ta.get_gene().get_matrix(k), 1) == 0);

#if 0
	int bin_count = 1;
	while(maj_indices[k].n_rows / (2 * bin_count) > mac_indices[k].n_rows)
	  bin_count++;

	std::cerr << "bin_count: " << bin_count << "\n";

	int nperm = ta.get_npermutations();

	mac_case_count[k] = ta.get_permute(k).random_case_count(
		nperm > 0 ? nperm : 1,
		mac_indices[k],
		maj_indices[k],
		ta.get_cov().get_probability(),
		ta.get_cov().get_ncases(),
		bin_count);

	for(const auto &v : mac_case_count[k]) {
	  maj_case_count[k].push_back(ta.get_cov().get_ncases() - v);
	}
#endif
	assert(mac_indices[k].n_rows + maj_indices[k].n_rows == ta.get_cov().get_nsamples());

	mac_odds[k] = ta.get_cov().get_odds()(mac_indices[k]);
	maj_odds[k] = ta.get_cov().get_odds()(maj_indices[k]);
  }

  int iter = 0;
  arma::vec phenotypes = ta.get_cov().get_phenotype_vector();

  if(ta.get_tp().permute_set) {
	pset_ofs.open(ta.get_tp().permute_set_path, std::ios_base::app);
  }

  while (iter < ta.get_npermutations()) {
	// For each transcript in the gene
	transcript_no = -1;
	for (auto &v : ta.results) {
	  std::vector<std::vector<int32_t>> permutations;
	  const std::string &k = v.second.transcript;

	  transcript_no++;

	  // SKAT corrects for covariates so we don't use this permutation approach
	  if (!(ta.get_methods().str() == "SKAT") && !(ta.get_methods().str() == "SKATO") && !(ta.get_methods().str() == "BURDEN")) {
		// Permute minor allele carriers
		// arma::vec temp_odds(mac_odds[k].n_rows + 1, arma::fill::zeros);
		// temp_odds(arma::span(0, mac_odds[k].n_rows - 1)) = mac_odds[k];
		// double group_odds = arma::mean(ta.get_cov().get_probability()(maj_indices[k])) / arma::mean(1 - ta.get_cov().get_probability()(maj_indices[k]));
		// temp_odds(temp_odds.n_rows - 1) = group_odds;

		// Get the number allotted to the large maj_allele_carrier bin
		// int nmaj_cases = permutations[0].back();
		// permutations[0].pop_back(); // Drop end
		// permutations[0].shrink_to_fit();

		permutations = ta.get_permute(k).permutations_maj_bin(1,
															  ta.get_cov().get_odds(),
															  ta.get_cov().get_ncases(),
															  mac_indices[k],
															  maj_indices[k]);
		arma::uword total_cases = 0;
		for (int i = 0; i < mac_indices[k].n_elem; i++) {
		  phenotypes(mac_indices[k](i)) = permutations[0][i];
		  if (permutations[0][i] == 1)
			total_cases++; // Count cases in mac
		}

		if(ta.get_tp().permute_set && transcript_no == 0) {
		  for(const auto &v : phenotypes) {
		    if(v > 1 || v < 0) {
		      throw(std::logic_error("phenotypes should be either 0 or 1"));
		    }
		    pset_ofs << v << "\t";
		  }
		  pset_ofs << "\n";
		}

		arma::vec temp(maj_indices[k].n_elem, arma::fill::zeros);
		arma::uword remaining_cases = ta.get_cov().get_ncases() - total_cases;
		temp(arma::span(0, remaining_cases - 1)) =
			arma::vec(ta.get_cov().get_ncases() - total_cases, arma::fill::ones);

		phenotypes(maj_indices[k]) = temp;
#if 0
		permutations = get_permutations(1, mac_odds[k], mac_case_count[k][iter]);

		phenotypes(mac_indices[k]) = arma::conv_to<arma::vec>::from(permutations[0]);

		arma::vec temp(maj_indices[k].n_rows, arma::fill::zeros);
		temp(arma::span(0, maj_case_count[k][iter] - 1)) = arma::vec(maj_case_count[k][iter], arma::fill::ones);
		phenotypes(maj_indices[k]) = arma::shuffle(temp);
#endif

#if 0
		if(maj_indices[k].n_rows > 500) {
		  arma::uword breaks = 10;
		  arma::uword split = maj_indices[k].n_rows / breaks;
		  std::vector<arma::span> span_vec;
		  std::vector<int32_t> bin_counts(breaks, 0);

		  // Sort odds to minimize variance in bins
		  arma::uvec odds_sort = arma::sort_index(maj_odds[k]);

		  // Get spans over breaks
		  for(arma::uword i = 0; i < breaks; i++) {
			if(i == breaks - 1) {
			  span_vec.emplace_back(arma::span(i * split, maj_indices[k].n_rows - 1));
			  bin_counts[i] = static_cast<int32_t>(maj_indices[k].n_rows - i * split);
			} else {
			  span_vec.emplace_back(arma::span(i * split, (i + 1) * split - 1));
			  bin_counts[i] = static_cast<int32_t>((i + 1) * split - i * split);
			}
		  }

		  arma::vec break_odds(breaks, arma::fill::zeros);
		  // Make temporaries of the sorted values to avoid templating problems
		  arma::vec odds_sorted = maj_odds[k](odds_sort);
		  arma::uvec indices_sorted = maj_indices[k](odds_sort);
		  // Get odds for group
		  for(arma::uword i = 0; i < breaks; i++) {
			break_odds(i) = arma::mean(odds_sorted(span_vec[i]));
		  }

		  // Get the number in each bin
		  std::vector<std::vector<int32_t>> span_cases = cases_in_bins(1, break_odds, maj_case_count[k][iter], bin_counts);

		  for(arma::uword i = 0; i < breaks; i++) {
			arma::vec temp_odds = odds_sorted(span_vec[i]);
			permutations = get_permutations(1, temp_odds, span_cases[0][i]);
			phenotypes(indices_sorted(span_vec[i])) = arma::conv_to<arma::vec>::from(permutations[0]);
		  }
		} else {
		  permutations = get_permutations(1, maj_odds[k], maj_case_count[k][iter]);
		  phenotypes(maj_indices[k]) = arma::conv_to<arma::vec>::from(permutations[0]);
		}

		// permutations = get_permutations(1, maj_odds[k], maj_case_count[k][iter]);
		// phenotypes(maj_indices[k]) = arma::conv_to<arma::vec>::from(permutations[0]);
		// for (arma::uword i = 0; i < maj_indices[k].n_rows; i++) {
		//   phenotypes(maj_indices[k](i)) = permutations[0][i];
		// }

#else

		// permutations = get_permutations(1, maj_odds[k], maj_case_count[k][iter]);
		// phenotypes(maj_indices[k]) = arma::conv_to<arma::vec>::from(permutations[0]);
		// Test alternative approach: simple shuffle of remaining cases
		// arma::vec temp(maj_indices[k].n_rows, arma::fill::zeros);
		// temp(arma::span(0, maj_case_count[k][iter] - 1)) = arma::vec(maj_case_count[k][iter], arma::fill::ones);
		// phenotypes(maj_indices[k]) = arma::shuffle(temp);

#endif
		ta.get_cov().set_phenotype_vector(phenotypes);

		perm_val = call_method(ta.get_methods(), ta.get_gene(), ta.get_cov(), k, ta.get_tp(), false, false);
	  } else {
		if (transcript_no == 0) {
		  perm_val = call_method(ta.get_methods(), ta.get_gene(), ta.get_cov(), k, ta.get_tp(), true, false);
		} else {
		  perm_val = call_method(ta.get_methods(), ta.get_gene(), ta.get_cov(), k, ta.get_tp(), false, false);
		}
	  }

	  // ta.increment_permuted(v.second.transcript, perm_val);
	  v.second.permuted.push_back(perm_val);

	  check_perm(ta.get_methods().str(), perm_val, ta.success_threshold, v);

	  // Update total number of permutations
	  v.second.permutations++;

	  // Track when we reached 30
	  if (v.second.successes == ta.success_threshold && v.second.min_success_at < 0) {
		v.second.min_success_at = v.second.permutations;
	  }
	}
	// Stop iterating if all transcripts are finished
	if (std::all_of(ta.results.cbegin(), ta.results.cend(), [&](const auto &v) { return v.second.done; }))
	  break;
	iter++;
  }
  if (verbose_) {
	for (const auto &v : ta.results) {
	  std::cerr << "Stage 2: " << ta.get_gene().get_gene() << "\t" << v.second.transcript << "\t";
	  std::cerr << v.second.original << std::endl;
	}
  }
  if (ta.get_tp().permute_set) {
    pset_ofs.close();
    std::exit(0);
  }

  for (auto &v : ta.results) {
	double empirical;
	double midp;
	if (v.second.min_success_at > 0 && v.second.successes == ta.success_threshold + 1) {
	  std::uniform_int_distribution<> dis(v.second.min_success_at, v.second.permutations);

	  empirical = v.second.successes / (1. + dis(gen_));
	  midp = v.second.mid_successes / (1 + dis(gen_));
	} else {
	  empirical = (1. + v.second.successes) / (1. + v.second.permutations);
	  midp = (1. + v.second.mid_successes) / (1. + v.second.permutations);
	}

	// Success on every iteration
	if (empirical > 1) {
	  empirical = 1;
	}

	if (midp > 1) {
	  midp = 1;
	}

	v.second.empirical_p = empirical;
	v.second.empirical_midp = midp;
  }
  ta.set_stage(Stage::Done);
  dispatch(ta);
}

size_t TaskQueue::get_nthreads() {
  return nthreads_;
}

void TaskQueue::duplicate() {
  results_detail_ = results_;
}

std::vector<TaskArgs> &TaskQueue::get_results_duplicate() {
  return results_detail_;
}

void TaskQueue::check_perm(const std::string &method,
						   double perm_val,
						   int success_threshold,
						   std::pair<const std::string, Result> &v) {
  // SKATO returns a pvalue so we need to reverse the successes
  if (method == "SKATO" || method == "SKAT") {
	if (perm_val <= v.second.original) {
	  if (v.second.successes < success_threshold) {
		v.second.successes++;
		if (perm_val == v.second.original) {
		  v.second.mid_successes += 0.5;
		} else {
		  v.second.mid_successes++;
		}
	  } else {
		v.second.successes++;
		if (perm_val == v.second.original) {
		  v.second.mid_successes += 0.5;
		} else {
		  v.second.mid_successes++;
		}
		v.second.done = true;
	  }
	}
  } else {
	if (perm_val >= v.second.original) {
	  if (v.second.successes < success_threshold) {
		v.second.successes++;
		if (perm_val == v.second.original) {
		  v.second.mid_successes += 0.5;
		} else {
		  v.second.mid_successes++;
		}
	  } else {
		v.second.successes++;
		if (perm_val == v.second.original) {
		  v.second.mid_successes += 0.5;
		} else {
		  v.second.mid_successes++;
		}
		v.second.done = true;
	  }
	}
  }
}

double TaskQueue::call_method(Methods &method, Gene &gene, Covariates &cov, const std::string &k, TaskParams &tp, bool shuffle, bool detail) {
  if(tp.method == "BURDEN") {
    return method.BURDEN(gene, k, shuffle, tp.a, tp.b);
  } else if(tp.method == "CALPHA") {
    return method.CALPHA(gene, cov, k);
  } else if(tp.method == "CMC") {
    return method.CMC(gene, cov, k, tp.maf);
  } else if(tp.method == "SKAT") {
    return method.SKATR(gene, k, shuffle, tp.a, tp.b, detail);
  } else if(tp.method == "SKATO") {
    return method.SKATRO(gene, k, shuffle, tp.a, tp.b, detail);
  } else if(tp.method == "VAAST") {
    return method.Vaast(gene, cov, k, tp.score_only_minor, tp.score_only_alternative, 2.0, tp.group_size, detail);
  } else if(tp.method == "VT") {
    return method.VT(gene, cov, k);
  } else if(tp.method == "WSS") {
    return method.WSS(gene, cov, k);
  }
}

