//
// Created by Bohlender,Ryan James on 10/2/18.
//

#include "vaast.hpp"

#include <algorithm>

VariantGroup::VariantGroup(arma::mat X,
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

  // arma::mat Xcollapse;
  if (group_threshold > 0) {
	arma::uvec score_sort = arma::stable_sort_index(variant_scores);
	X = X.cols(score_sort);
	weights = weights(score_sort);

#if 0
	// Do grouping
	arma::uvec rare = arma::find(arma::sum(X) <= group_threshold);
	arma::uvec common = arma::find(arma::sum(X) > group_threshold);

	for (arma::uword i = 1; i < rare.n_rows; i++) {
	  arma::span tail(rare.n_rows - 1 - i, rare.n_rows - 1);
	  arma::span head(0, rare.n_rows - 1 - i);

	  arma::uvec singles = arma::join_vert(common, rare(head));
	  arma::uvec collapsed = rare(tail);

	  if (i == rare.n_rows - 1) {
		Xcollapse = arma::sum(X.cols(collapsed), 1);
		weight = arma::prod(weights(collapsed), 0);
	  } else {
		Xcollapse = arma::join_horiz(X.cols(singles), arma::sum(X.cols(collapsed), 1));
		weight = arma::join_vert(weights(singles), arma::prod(weights(collapsed), 0));
	  }

	  case_allele1 = Xcollapse.t() * Y;
	  control_allele1 = Xcollapse.t() * (1 - Y);
	  case_allele0 = 2 * n_case - case_allele1;
	  control_allele0 = 2 * n_control - control_allele1;

	  double new_score = Score();
	  if (new_score > score) {
		score = new_score;
	  } else {
		// Reset to best group
		i--;
		tail = arma::span(rare.n_rows - 1 - i, rare.n_rows - 1);
		head = arma::span(0, rare.n_rows - 1 - i);

		singles = arma::join_vert(common, rare(head));
		collapsed = rare(tail);

		if (i == rare.n_rows - 1) {
		  Xcollapse = arma::sum(X.cols(collapsed), 1);
		  weight = arma::prod(weights(collapsed), 0);
		} else {
		  Xcollapse = arma::join_horiz(X.cols(singles), arma::sum(X.cols(collapsed), 1));
		  weight = arma::join_vert(weights(singles), arma::prod(weights(collapsed), 0));
		}

		case_allele1 = Xcollapse.t() * Y;
		control_allele1 = Xcollapse.t() * (1 - Y);
		case_allele0 = 2 * n_case - case_allele1;
		control_allele0 = 2 * n_control - control_allele1;

		break;
	  }
	  nvariants = Xcollapse.n_cols; // Final number of variants
	  score = new_score > score ? new_score : score;
	}
#endif
	double last_score = -1;
	arma::sword i = X.n_cols - 1;
	arma::span merged;
	for (; i >= 0; i--) {
	  merged = arma::span(i, X.n_cols - 1);
	  arma::mat Xnew = arma::sum(X.cols(merged), 1);
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
	if (i == -1) {
	  // All variants merged
	  Xcollapse = arma::sum(X, 1);
	  weight = arma::prod(weights, 0);

	  score = Score(Xcollapse, Y, weight);

	} else {
	  merged = arma::span(i + 1, X.n_cols - 1);
	  arma::span unmerged(0, i);

	  Xcollapse = arma::join_horiz(X.cols(unmerged), arma::sum(X.cols(merged), 1));
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

double VariantGroup::Score(const arma::mat &X, const arma::vec &Y, const arma::vec &w) {
  case_allele1 = X.t() * Y;
  control_allele1 = X.t() * (1 - Y);
  case_allele0 = 2 * n_case - case_allele1;
  control_allele0 = 2 * n_control - control_allele1;

  nvariants = X.n_cols; // Number of variants

  arma::vec log_lh = LRT();
  variant_scores = 2.0 * (log_lh + arma::log(w)) - site_penalty;

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
			 Covariates &cov,
			 const std::string &k,
			 bool score_only_minor,
			 bool score_only_alternative,
			 double site_penalty,
			 arma::uword group_threshold,
			 bool detail)
	: som(score_only_minor), soa(score_only_alternative), detail(detail), k(k), group_threshold(group_threshold),
	  site_penalty(site_penalty),
	  X(gene.get_matrix(k)), Y(cov.get_phenotype_vector()) {
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
		for (arma::uword i = 0; i < X.n_cols; i++) {
		  arma::mat Xnew = X;
		  arma::vec Wnew = weights;
		  std::vector<std::string> positions = gene.get_positions(k);

		  positions.erase(positions.begin() + i);
		  Xnew.shed_col(i);
		  Wnew.shed_row(i);

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

double VAAST::Score() {
  double val = 0;
  for (const auto &g : groups) {
    val += arma::accu(g.variant_scores);
  }
  return val;
}

double VAAST::Score(const arma::mat &X, const arma::vec &Y, const arma::vec &w) {
  // Get sample sizes
  n_case = static_cast<arma::uword>(arma::sum(Y));
  n_control = static_cast<arma::uword>(arma::sum(1 - Y));

  // Get carrier counts
  case_allele1 = X.t() * Y;
  control_allele1 = X.t() * (1 - Y);
  case_allele0 = 2. * n_case - case_allele1;
  control_allele0 = 2. * n_control - control_allele1;

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

void VAAST::variant_grouping(const arma::mat &X,
							 const arma::vec &Y,
							 const arma::vec &w,
							 std::vector<std::string> &positions) {
  groups.clear();
  if (group_threshold == 0) {
	groups.push_back(VariantGroup(X, Y, w, 0, site_penalty, som, soa));
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
	arma::mat Xnew = X.cols(idx);
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
	  groups.push_back(VariantGroup(Xnew, Y, Wvec, 0, site_penalty, som, soa));
	  continue;
	}

	arma::mat Xcollapsible = Xnew.cols(target);
	arma::vec Wcollapsible = Wvec(target);
	arma::vec score_set = vaast_site_scores(idx).eval()(target);
	arma::uvec score_sort = arma::sort_index(score_set);

	arma::uvec weight_sort = arma::sort_index(Wcollapsible);

	arma::uword ngroups;
	if (type.first == "SNV") {
	  groups.push_back(VariantGroup(Xcollapsible.cols(score_sort), Y, Wcollapsible(score_sort), group_threshold, site_penalty, som, soa));
	  if (remaining.n_elem > 0) {
		groups.push_back(VariantGroup(Xnew.cols(remaining), Y, Wvec(remaining), 0, site_penalty, som, soa));
	  }
	} else {
	  ngroups = static_cast<arma::uword>(std::ceil(weight_sort.n_elem / 5.));
	  for (arma::uword i = 0; i < ngroups; i++) {
	    // Sort on weight for grouping
		arma::span cur_span(i * 5, std::min(i * 5 + 4, weight_sort.n_elem - 1));
		arma::mat Xmat = Xcollapsible.cols(weight_sort(cur_span));
		arma::vec weight_spanned = Wcollapsible(weight_sort(cur_span));
		// Resort by score for merging
		arma::vec score_spanned = score_set(weight_sort(cur_span));
		score_sort = arma::sort_index(score_spanned);
		groups.push_back(VariantGroup(Xmat.cols(score_sort), Y, weight_spanned(score_sort), group_threshold, site_penalty, som, soa));
	  }
	  // Don't group the remaining variants.
	  if (remaining.n_elem > 0) {
		groups.push_back(VariantGroup(Xnew.cols(remaining), Y, Wvec(remaining), 0, site_penalty, som, soa));
	  }
	}
  }
}

//! \brief Getter for the VAAST score
//! \return double containing the VAAST score
double VAAST::get_score() {
  return score;
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

void VAAST::variant_bitmask(const arma::mat &X, const arma::vec &Y, const arma::vec &w) {
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
