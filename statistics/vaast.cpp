//
// Created by Bohlender,Ryan James on 10/2/18.
//

#include "vaast.hpp"

#include <algorithm>

VariantGroup::VariantGroup(arma::sp_mat X,
						   arma::vec Y,
						   arma::vec weights,
						   arma::uword group_threshold,
						   double site_penalty,
						   bool score_only_minor,
						   bool score_only_alternative)
	: Xcollapse(X), Y(Y), site_penalty(site_penalty), som(score_only_minor), soa(score_only_alternative) {
  n_case = static_cast<arma::uword>(arma::sum(Y));
  n_control = static_cast<arma::uword>(arma::sum(1 - Y));

  // Get initial score:
  weight = weights;

  score = Score(X, Y, weights);

  arma::uvec score_sort = arma::stable_sort_index(variant_scores);
  arma::sp_mat Xmat(X.n_rows, X.n_cols);
  // X = X.cols(score_sort);
  arma::uword i = 0;
  for(const auto &v : score_sort) {
	Xmat(arma::span::all, i) = X.col(v);
	i++;
  }

  // arma::mat Xcollapse;
  if (group_threshold > 0) {
	weights = weights(score_sort);

	double last_score = -1;
	arma::sword j = X.n_cols - 1;
	arma::span merged;
	for (; j >= 0; j--) {
	  merged = arma::span(j, X.n_cols - 1);
	  arma::sp_mat Xnew = arma::sum(Xmat.cols(j, X.n_cols - 1), 1);
	  arma::vec Wnew = arma::prod(weights(merged), 0);

	  double new_score = Score(Xnew, Y, Wnew);
	  if (new_score > last_score || last_score < 0) {
		last_score = new_score;
		continue;
	  } else {
		break;
	  }
	}
	// Finish grouping
	if (j == -1) {
	  // All variants merged
	  Xcollapse = arma::sum(Xmat, 1);
	  weight = arma::prod(weights, 0);

	  score = Score(Xcollapse, Y, weight);

	} else {
	  merged = arma::span(j + 1, X.n_cols - 1);
	  arma::span unmerged(0, j);

	  Xcollapse = arma::join_horiz(Xmat.cols(0, j), arma::sum(X.cols(j + 1, X.n_cols - 1), 1));
	  weight = arma::join_vert(weights(unmerged), arma::prod(weights, 0));

	  score = Score(Xcollapse, Y, weight);
	}
  } else {
	// Keep initial score
	return;
  }

}

void VariantGroup::variant_mask() {
#if 1
  arma::uvec tmpmask = variant_scores <= 0;
  arma::uvec scenario1 =
	  ((control_allele1 > control_allele0) + (case_allele1 % control_allele0 < case_allele0 % control_allele1)) >= 1;
  mask = arma::find(tmpmask + scenario1 >= 1);
#else
  arma::uvec tmpmask = arma::find(variant_scores <= 0);
  std::vector<arma::uword> scenario;
  // mask sites where major allele is more common in cases
  /*
	# scenario 1: 1 is minor allele and it's more frequent in case
		scenario1= (    (control_allele1 <= control_allele0 )
			&   (case_allele1 * control_allele0 >= case_allele0 *
				control_allele1)
		)
	# scenario 2: 0 is minor allele and it's more frequent in case
		scenario2= (    (control_allele1 >= control_allele0 )
			&   (case_allele1 * control_allele0 <= case_allele0 *
				control_allele1)
		)
	*/
  case_allele1 = Xcollapse.t() * Y;
  control_allele1 = Xcollapse.t() * (1 - Y);
  case_allele0 = 2. * n_case - case_allele1;
  control_allele0 = 2. * n_control - control_allele1;

  for (arma::uword i = 0; i < case_allele1.n_elem; i++) {
	bool in_mask = false;
	for (arma::uword j = 0; j < tmpmask.size(); j++) {
	  if (i == tmpmask(j)) {
		in_mask = true;
		break;
	  }
	}
	// Merge
	if (in_mask) {
	  scenario.push_back(i);
	  continue;
	}
	if (som) {
	  // If both scenarios are false, append to mask
	  bool scenario1 = !((control_allele1(i) <= control_allele0(i))
		  & ((case_allele1(i) * control_allele0(i) >= case_allele0(i) * control_allele1(i))));
	  bool scenario2 = !((control_allele1(i) >= control_allele0(i))
		  & ((case_allele1(i) * control_allele0(i) <= case_allele0(i) * control_allele1(i))));

	  if (soa) {
		if (scenario1)
		  scenario.push_back(i);
	  } else {
		if (scenario1 & scenario2)
		  scenario.push_back(i);
	  }
	}
  }
  mask = arma::conv_to<arma::uvec>::from(scenario);
#endif
}

double VariantGroup::Score(const arma::sp_mat &X, const arma::vec &Y, const arma::vec &w) {
  case_allele1 = X.t() * Y;
  control_allele1 = X.t() * (1 - Y);
  case_allele0 = 2 * n_case - case_allele1;
  control_allele0 = 2 * n_control - control_allele1;

  nvariants = X.n_cols; // Number of variants

  arma::vec log_lh = LRT();
  variant_scores = 2.0 * (log_lh + arma::log(w));
  for(auto & v : variant_scores) {
    if (v > site_penalty)
      v -= site_penalty;
  }

  variant_mask();
  variant_scores(mask).zeros();

  double val = arma::accu(variant_scores);
  return (val >= 0) ? val : 0; // Mask all negative values
}

arma::vec VariantGroup::LRT() {
  arma::vec alt_control_freq = control_allele1 / (control_allele0 + control_allele1);
  arma::vec alt_case_freq = case_allele1 / (case_allele0 + case_allele1);

  arma::vec
	  null_freq = (case_allele1 + control_allele1) / (case_allele0 + case_allele1 + control_allele0 + control_allele1);

  arma::vec alt_log_lh = log_likelihood(alt_control_freq, control_allele0, control_allele1)
	  + log_likelihood(alt_case_freq, case_allele0, case_allele1);
  arma::vec null_log_lh = log_likelihood(null_freq, control_allele0, control_allele1)
	  + log_likelihood(null_freq, case_allele0, case_allele1);

  return alt_log_lh - null_log_lh;
}

arma::vec VariantGroup::log_likelihood(arma::vec &freq, arma::vec &allele0, arma::vec &allele1) {
  // Prevent numerical issues
  arma::vec clamped = arma::clamp(freq, 1e-9, 1.0 - 1e-9);

  return allele1 % arma::log(clamped) + allele0 % arma::log(1.0 - clamped);
}

VAAST::VAAST(Gene &gene,
			 arma::vec &Y,
			 const std::string &k,
			 bool score_only_minor,
			 bool score_only_alternative,
			 double site_penalty,
			 arma::uword group_threshold,
			 bool detail,
			 bool biallelic)
	: som(score_only_minor), soa(score_only_alternative), detail(detail), biallelic(biallelic), k(k), group_threshold(group_threshold),
	  site_penalty(site_penalty),
	  X(gene.get_matrix(k)), Y(Y) {
  // Verify weights are okay
  check_weights(gene);

  if (group_threshold == 0 || X.n_cols == 1) {
	score = Score(X, Y, weights);
	if (detail)
	  gene.set_scores(k, vaast_site_scores);
  } else {
	score = Score(X, Y, weights);
	variant_grouping(X, Y, weights, gene.get_positions(k));
	score = Score();

  }

  if (detail) {
	if (group_threshold == 0) {
	  gene.set_scores(k, vaast_site_scores);
	} else {
	  expanded_scores.reshape(X.n_cols, 1);
	  if (X.n_cols > 1) {
		arma::sp_mat Xnew(X);
		arma::vec Wnew(weights);
		std::vector<std::string> positions = gene.get_positions(k);
		for (arma::uword i = 0; i < X.n_cols; i++) {
		  if(i == 0) {
			positions.erase(positions.begin());
			Xnew.shed_col(i);
			Wnew.shed_row(i);
		  } else {
			positions[i - 1] = gene.get_positions(k)[i - 1];
			Xnew.col(i - 1) = X.col(i - 1);
			Wnew(i - 1) = weights(i - 1);
		  }

		  Score(Xnew, Y, Wnew); // Do normal individual variant scoring
		  variant_grouping(Xnew, Y, Wnew, positions); // Redo grouping
		  double new_score = Score();
		  expanded_scores(i) = (score - new_score >= 0) ? score - new_score : 0; // Calculate score
		}
		gene.set_scores(k, expanded_scores);
	  } else {
		expanded_scores(0) = vaast_site_scores(0);
		gene.set_scores(k, expanded_scores);
	  }
	}
  }
}

VAAST::VAAST(arma::sp_mat X,
			 arma::vec &Y,
			 arma::vec &weights,
			 std::vector<std::string> &positions_,
			 const std::string &k,
			 bool score_only_minor,
			 bool score_only_alternative,
			 bool biallelic,
			 arma::uword group_threshold,
			 double site_penalty)
	: som(score_only_minor), soa(score_only_alternative), detail(true), biallelic(biallelic), k(k), group_threshold(group_threshold),
	  site_penalty(site_penalty),
	  X(std::move(X)), Y(Y) {

  if (group_threshold == 0 || X.n_cols == 1) {
	score = Score(X, Y, weights);
  } else {
	score = Score(X, Y, weights);
	variant_grouping(X, Y, weights, positions_);
	score = Score();

  }

  if(group_threshold > 0) {
	expanded_scores.reshape(X.n_cols, 1);
	if (X.n_cols > 1) {
	  arma::sp_mat Xnew(X);
	  arma::vec Wnew(weights);
	  std::vector<std::string> positions = positions_;
	  for (arma::uword i = 0; i < X.n_cols; i++) {
	    if(i == 0) {
		  positions.erase(positions.begin());
		  Xnew.shed_col(i);
		  Wnew.shed_row(i);
	    } else {
	      positions[i - 1] = positions_[i - 1];
	      Xnew.col(i - 1) = X.col(i - 1);
	      Wnew(i - 1) = weights(i - 1);
	    }
#if 0
		for(auto it = X.begin(); it != X.end(); it++) {
		  if(it.col() == i) {
		    continue;
		  } else if(it.col() < i && *it > 0) {
		    Xnew(it.row(), it.col()) = *it;
		  } else if(it.col() > i && *it > 0) {
		    Xnew(it.row(), it.col() - 1) = *it;
		  }
		}
		for(auto it = weights.begin(); it != weights.end(); it++) {
		  auto j = std::distance(weights.begin(), it);
		  if(i == j) {
		    continue;
		  } else if(j < i) {
		    Wnew(j) = *it;
		  } else if(j > i) {
		    Wnew(j - 1) = *it;
		  }
		}
#endif

		Score(Xnew, Y, Wnew); // Do normal individual variant scoring
		variant_grouping(Xnew, Y, Wnew, positions); // Redo grouping
		double new_score = Score();
		expanded_scores(i) = (score - new_score >= 0) ? score - new_score : 0; // Calculate score
	  }
	} else {
	  expanded_scores(0) = vaast_site_scores(0);
	}
  } else {
    expanded_scores = vaast_site_scores;
  }
}

double VAAST::Score() {
  double val = 0;
  for (const auto &g : groups) {
    val += arma::accu(g.variant_scores);
  }
  return val;
}

double VAAST::Score(const arma::sp_mat &X, const arma::vec &Y, const arma::vec &w) {
  // Get sample sizes
  n_case = static_cast<arma::uword>(arma::sum(Y));
  n_control = static_cast<arma::uword>(arma::sum(1 - Y));

  // Get carrier counts
  case_allele1 = X.t() * Y;
  control_allele1 = X.t() * (1 - Y);
  case_allele0 = 2. * n_case - case_allele1;
  control_allele0 = 2. * n_control - control_allele1;

  // Handle collapsed variant
  if(biallelic) {
    // Via email:
    // The idea is to add an additional term to the VAAST CLRT that captures the frequency of biallelic variants in cases and controls.
    // For starters, we can try collapsing all variants with a positive VAAST score into one biallelic term.
	arma::vec log_lh = LRT();
	vaast_site_scores = 2.0 * (log_lh + arma::log(w)) - site_penalty;
	variant_bitmask(X, Y, w);
	vaast_site_scores(mask).zeros();

	// Get all variants with positive vaast score
	arma::uvec collapse = arma::find(vaast_site_scores > 0);
	arma::vec variant(n_case + n_control, arma::fill::zeros);

	for(arma::uword i = 0; i < collapse.n_elem; i++) {
	  for(arma::uword j = 0; j < n_case + n_control; j++) {
	    if(X(j, i) == 1) {
	      variant(j) += 1;
	    }
	    if(X(j, i) == 2) {
		  variant(j) += 2;
	    }
	  }
	}
	variant(arma::find(variant == 1)).fill(0); // Set all single hets to 0
	variant(arma::find(variant > 1)).fill(1); // Set all compound hets and homozygous rare carriers to 1

	case_allele1 = variant.t() * Y;
	control_allele1 = variant.t() * (1 - Y);
	case_allele0 = n_case - case_allele1;
	control_allele0 = n_control - control_allele1; // Not 2 * because we're only allowing presence absence.

	arma::vec biallelic_term = LRT();
	double val = arma::accu(vaast_site_scores) + arma::accu(biallelic_term);
	return (val >= 0) ? val : 0; // Mask all negative values
  }

  arma::vec log_lh = LRT();
  vaast_site_scores = 2.0 * (log_lh + arma::log(w)) - site_penalty;
  variant_bitmask(X, Y, w);
  vaast_site_scores(mask).zeros();

  double val = arma::accu(vaast_site_scores);
  return (val >= 0) ? val : 0; // Mask all negative values
}

arma::vec VAAST::LRT() {
  arma::vec alt_control_freq = control_allele1 / (control_allele0 + control_allele1);
  arma::vec alt_case_freq = case_allele1 / (case_allele0 + case_allele1);

  arma::vec
	  null_freq = (case_allele1 + control_allele1) / (case_allele0 + case_allele1 + control_allele0 + control_allele1);

  arma::vec alt_log_lh = log_likelihood(alt_control_freq, control_allele0, control_allele1)
	  + log_likelihood(alt_case_freq, case_allele0, case_allele1);
  arma::vec null_log_lh = log_likelihood(null_freq, control_allele0, control_allele1)
	  + log_likelihood(null_freq, case_allele0, case_allele1);

  return alt_log_lh - null_log_lh;
}

arma::vec VAAST::log_likelihood(arma::vec &freq, arma::vec &allele0, arma::vec &allele1) {
  // Prevent numerical issues
  arma::vec clamped = arma::clamp(freq, 1e-9, 1.0 - 1e-9);

  return allele1 % arma::log(clamped) + allele0 % arma::log(1.0 - clamped);
}

void VAAST::check_weights(Gene &gene) {
  if (gene.is_weighted(k)) {
	weights = gene.get_weights(k);
	return;
  }

  weights = arma::vec(X.n_cols, arma::fill::ones);
  gene.set_weights(k, weights);
}

void VAAST::variant_grouping(const arma::sp_mat &X,
							 const arma::vec &Y,
							 const arma::vec &w,
							 std::vector<std::string> &positions) {
  groups.clear();
  if (group_threshold == 0) {
	groups.emplace_back(VariantGroup(X, Y, w, 0, site_penalty, som, soa));
	return;
  }
  // Initially group into indels vs. snvs
  std::vector<arma::uword> insertion;
  std::vector<arma::uword> deletion;
  std::vector<arma::uword> mnp;
  std::vector<arma::uword> SNVs;
  arma::uword i = 0;
  for (const auto &p : positions) {
	RJBUtil::Splitter<std::string> splitter(p, "-");
	if (splitter.back() == "SNV") {
	  SNVs.push_back(i);
	} else if (splitter.back() == "insertion") {
	  insertion.push_back(i);
	} else if (splitter.back() == "deletion") {
	  deletion.push_back(i);
	} else if (splitter.back() == "complex_substitution") {
	  mnp.push_back(i);
	} else {
	  throw(std::logic_error("Wrong variant type."));
	}
	i++;
  }
  // Identify rare variants <= group_threshold
  // More common in cases than controls
  // All other variants are singletons
  std::map<std::string, arma::uvec> type_idx_map;
  type_idx_map["insertion"] = arma::conv_to<arma::uvec>::from(insertion);
  type_idx_map["deletion"] = arma::conv_to<arma::uvec>::from(deletion);
  type_idx_map["mnp"] = arma::conv_to<arma::uvec>::from(mnp);
  type_idx_map["SNV"] = arma::conv_to<arma::uvec>::from(SNVs);

  for (const auto &type : type_idx_map) {
	arma::uvec idx = type.second;
	if (idx.n_elem == 0) {
	  continue;
	}
	// arma::sp_mat Xnew = X.cols(idx);

	arma::sp_mat Xnew(X.n_rows, idx.n_elem);
	i = 0;
	for(const auto &v : idx) {
	  Xnew.col(i) = X.col(v);
	  i++;
	}
	arma::vec Wvec = w(idx);
	arma::uvec indices = arma::regspace<arma::uvec>(0, Xnew.n_cols - 1);

	arma::vec case_count1 = Xnew.t() * Y;
	arma::vec cont_count1 = Xnew.t() * (1 - Y);

	arma::uvec rare = arma::find((case_count1 <= group_threshold) && (vaast_site_scores(idx) > 0.));

	// Variants more common in cases than in remaining.
	arma::uvec case_frequent = arma::find(case_count1 > cont_count1);

	arma::uvec target = arma::intersect(rare, case_frequent);
	arma::uvec remaining = setdiff(indices, target);

	// Nothing to group; finish early
	if (target.n_elem == 0) {
	  groups.emplace_back(VariantGroup(Xnew, Y, Wvec, 0, site_penalty, som, soa));
	  continue;
	}

	arma::sp_mat Xcollapsible(Xnew.n_rows, target.n_elem);
	i = 0;
	for(const auto &v : target) {
	  Xcollapsible.col(i) = Xnew.col(v);
	  i++;
	}
	arma::vec Wcollapsible = Wvec(target);
	arma::vec score_set = vaast_site_scores(idx).eval()(target);
	arma::uvec score_sort = arma::sort_index(score_set);

	arma::uvec weight_sort = arma::sort_index(Wcollapsible);

	arma::uword ngroups;
	if (type.first == "SNV") {
	  groups.emplace_back(VariantGroup(Xcollapsible, Y, Wcollapsible(score_sort), group_threshold, site_penalty, som, soa));
	  if (remaining.n_elem > 0) {
	    arma::sp_mat Xremain(Xnew.n_rows, remaining.n_elem);
	    i = 0;
	    for(const auto &v : remaining) {
	      Xremain.col(i) = Xnew.col(v);
		  i++;
	    }
		groups.emplace_back(VariantGroup(Xremain, Y, Wvec(remaining), 0, site_penalty, som, soa));
	  }
	} else {
	  ngroups = static_cast<arma::uword>(std::ceil(weight_sort.n_elem / 5.));
	  for (i = 0; i < ngroups; i++) {
	    // Sort on weight for grouping
		arma::span cur_span(i * 5, std::min(i * 5 + 4, weight_sort.n_elem - 1));
		arma::sp_mat Xmat(Xcollapsible.n_rows, cur_span.b - cur_span.a + 1);
		arma::uword j = 0;
		for(const auto &v : weight_sort(cur_span)) {
		  Xmat.col(j) = Xcollapsible.col(v);
		  j++;
		}
		arma::vec weight_spanned = Wcollapsible(weight_sort(cur_span));
		// Resort by score for merging
		arma::vec score_spanned = score_set(weight_sort(cur_span));
		score_sort = arma::sort_index(score_spanned);
		groups.emplace_back(VariantGroup(Xmat, Y, weight_spanned(score_sort), group_threshold, site_penalty, som, soa));
	  }
	  // Don't group the remaining variants.
	  if (remaining.n_elem > 0) {
		arma::sp_mat Xremain(Xnew.n_rows, remaining.n_elem);
		i = 0;
		for(const auto &v : remaining) {
		  Xremain.col(i) = Xnew.col(v);
		  i++;
		}
		groups.emplace_back(VariantGroup(Xremain, Y, Wvec(remaining), 0, site_penalty, som, soa));
	  }
	}
  }
}

//! \brief Calculate the set difference between two uvec
//! \param x uvec with elements to keep
//! \param y uvec with elements to exclude
//! \return copy of uvec excluding any element in y
arma::uvec VAAST::setdiff(arma::uvec x, arma::uvec y) {
  for (size_t j = 0; j < y.n_elem; j++) {
	arma::uvec q1 = arma::find(x == y[j]);
	if (!q1.empty()) {
	  x.shed_row(q1(0));
	}
  }

  return x;
}

void VAAST::variant_bitmask(const arma::sp_mat &X, const arma::vec &Y, const arma::vec &w) {
  // mask sites where major allele is more common in cases
  // We can simplify because the 1 is always the minor allele
  // arma::uvec tmpmask = vaast_site_scores <= 0;
  // arma::umat scenario1 = case_allele1 > control_allele1;
  // arma::umat scenario2 =

  /*
	# scenario 1: 1 is minor allele and it's more frequent in case
		scenario1= (    (control_allele1 <= control_allele0 )
			&   (case_allele1 * control_allele0 >= case_allele0 *
				control_allele1)
		)
	# scenario 2: 0 is minor allele and it's more frequent in case
		scenario2= (    (control_allele1 >= control_allele0 )
			&   (case_allele1 * control_allele0 <= case_allele0 *
				control_allele1)
		)
    */
#if 1
  arma::uvec tmpmask = vaast_site_scores <= 0;
  arma::uvec scenario1 =
	  ((control_allele1 > control_allele0) + (case_allele1 % control_allele0 < case_allele0 % control_allele1)) >= 1;
  mask = arma::find(tmpmask + scenario1 >= 1);
#else
  arma::uvec tmpmask = arma::find(vaast_site_scores <= 0);
  std::vector<arma::uword> scenario;
  for (arma::uword i = 0; i < case_allele1.n_elem; i++) {
	bool in_mask = false;
	for (arma::uword j = 0; j < tmpmask.size(); j++) {
	  if (i == tmpmask(j)) {
		in_mask = true;
		break;
	  }
	}
	// Merge
	if (in_mask) {
	  scenario.push_back(i);
	  continue;
	}
	if (som) {
	  // If both scenarios are false, append to mask
	  bool scenario1 = !((control_allele1(i) <= control_allele0(i))
		  & ((case_allele1(i) * control_allele0(i) >= case_allele0(i) * control_allele1(i))));
	  bool scenario2 = !((control_allele1(i) >= control_allele0(i))
		  & ((case_allele1(i) * control_allele0(i) <= case_allele0(i) * control_allele1(i))));

	  if (soa) {
		if (scenario1)
		  scenario.push_back(i);
	  } else {
		if (scenario1 & scenario2)
		  scenario.push_back(i);
	  }
	}
  }
  mask = arma::conv_to<arma::uvec>::from(scenario);
#endif
}
