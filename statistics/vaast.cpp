//
// Created by Bohlender,Ryan James on 10/2/18.
//

#include "vaast.hpp"

#include <algorithm>

VariantGroup::VariantGroup(arma::mat X, arma::vec Y, arma::vec weights, arma::uword group_threshold, double site_penalty)
: X(X), Y(Y), site_penalty(site_penalty) {
  double n_case = arma::sum(Y);
  double n_control = arma::sum(1 - Y);

  // Get initial score:
  weight = weights;

  case_allele1 = X.t() * Y;
  control_allele1 = X.t() * (1 - Y);
  case_allele0 = 2 * n_case - case_allele1;
  control_allele0 = 2 * n_control - control_allele1;

  score = Score();


  // Do grouping
  arma::uvec rare = arma::find(arma::sum(X) <= group_threshold);
  arma::uvec common = arma::find(arma::sum(X) > group_threshold);

  arma::mat Xcollapse;
  if(rare.n_rows >= 2 || group_threshold > 0) {
    for(arma::uword i = 1; i < rare.n_rows; i++) {
      arma::span tail(rare.n_rows - 1 - i, rare.n_rows - 1);
      arma::span head(0, rare.n_rows - 1 - i);

      arma::uvec singles = arma::join_vert(common, rare(head));
      arma::uvec collapsed = rare(tail);

      if(i == rare.n_rows - 1) {
		Xcollapse = arma::sum(X.cols(collapsed), 1);
		weight = arma::sum(weights(collapsed), 0);
      } else {
		Xcollapse = arma::join_horiz(X.cols(singles), arma::sum(X.cols(collapsed), 1));
		weight = arma::join_vert(weights(singles), arma::sum(weights(collapsed), 0));
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

		if(i == rare.n_rows - 1) {
		  Xcollapse = arma::sum(X.cols(collapsed), 1);
		  weight = arma::sum(weights(collapsed), 0);
		} else {
		  Xcollapse = arma::join_horiz(X.cols(singles), arma::sum(X.cols(collapsed), 1));
		  weight = arma::join_vert(weights(singles), arma::sum(weights(collapsed), 0));
		}

		case_allele1 = Xcollapse.t() * Y;
		control_allele1 = Xcollapse.t() * (1 - Y);
		case_allele0 = 2 * n_case - case_allele1;
		control_allele0 = 2 * n_control - control_allele1;
	    break;
	  }
	  score = new_score > score ? new_score : score;
    }
  } else {
    // Keep initial score
    return;
  }

}

void VariantGroup::mask() {
  score = 0;
}

double VariantGroup::Score() {
  arma::vec log_lh = LRT();
  variant_scores = 2.0 * (log_lh + arma::log(weight)) - site_penalty;
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

#if 0
double VariantGroup::LRT() {
  double alt_control_freq = control_allele1 / (control_allele0 + control_allele1);
  double alt_case_freq = case_allele1 / (case_allele0 + case_allele1);

  double null_freq = (case_allele1 + control_allele1) / (case_allele0 + case_allele1 + control_allele0 + control_allele1);

  double alt_log_lh = log_likelihood(alt_control_freq, control_allele0, control_allele1)
	  + log_likelihood(alt_case_freq, case_allele0, case_allele1);
  double null_log_lh = log_likelihood(null_freq, control_allele0, control_allele1)
	  + log_likelihood(null_freq, case_allele0, case_allele1);

  return alt_log_lh - null_log_lh;
}

double VariantGroup::log_likelihood(double freq, double allele0, double allele1) {
  // Prevent numerical issues
  auto clamp = [](double x, double lower, double upper) {
	if (x > upper) {
	  return upper;
	} else if(x < lower) {
	  return lower;
	} else {
	  return x;
	}
  };

  double clamped = clamp(freq, 1e-9, 1.0 - 1e-9);

  return allele1 * std::log(clamped) + allele0 * std::log(1.0 - clamped);
}
#endif

VAAST::VAAST(Gene &gene,
			 Covariates &cov,
			 const std::string &k,
			 bool score_only_minor,
			 bool score_only_alternative,
			 double site_penalty,
			 arma::uword group_threshold,
			 bool detail)
	: som(score_only_minor), soa(score_only_alternative), detail(detail), k(k), group_threshold(group_threshold), site_penalty(site_penalty),
	  X(gene.get_matrix(k)), Y(cov.get_phenotype_vector()) {
  // Verify weights are okay
  check_weights(gene);

  // Get sample sizes
  n_case = static_cast<arma::uword>(arma::sum(Y));
  n_control = static_cast<arma::uword>(arma::sum(1 - Y));

  // Get carrier counts
  case_allele1 = X.t() * Y;
  control_allele1 = X.t() * (1 - Y);
  case_allele0 = 2. * n_case - case_allele1;
  control_allele0 = 2. * n_control - control_allele1;

  variant_grouping(X, Y, weights);
  score = Score(X, Y, weights);

  if(detail) {
	expanded_scores.reshape(X.n_cols, 1);
	for(arma::uword i = 0; i < X.n_cols; i++) {
	  arma::mat Xnew = X;
	  arma::vec Wnew = weights;

	  Xnew.shed_col(i);
	  Wnew.shed_row(i);

	  variant_grouping(Xnew, Y, Wnew); // Redo grouping
	  double new_score = Score(Xnew, Y, Wnew);
	  expanded_scores(i) = (score - new_score >= 0) ? score - new_score : 0; // Calculate score
	}
	gene.set_scores(k, expanded_scores);
  }
}

double VAAST::Score(const arma::mat &X, const arma::vec &Y, const arma::vec &w) {
  arma::uword i = 0;
  vaast_site_scores.reshape(groups.size(), 1);
  for(auto &v : groups) {
	vaast_site_scores(i) = v.score;
	i++;
  }
  variant_bitmask(X, Y, w); // Update mask
  for(auto &i : mask) {
	groups[i].mask();
  }
  return arma::sum(vaast_site_scores);
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

void VAAST::variant_grouping(const arma::mat &X, const arma::vec &Y, const arma::vec &w) {
  groups.clear();
  // Identify rare variants <= group_threshold
  // More common in cases than controls
  // All other variants are singletons
  arma::uvec indices = arma::regspace<arma::uvec>(0, X.n_cols - 1);

  arma::uvec rare = arma::find(arma::sum(X) <= group_threshold);
  arma::uvec common = arma::find(arma::sum(X) > group_threshold);

  arma::vec case_count1 = X.cols(rare).t() * Y;
  arma::vec cont_count1 = X.cols(rare).t() * (1 - Y);

  arma::uvec case_frequent = arma::find(case_count1 > cont_count1);

  arma::uvec target = arma::intersect(rare, case_frequent);
  arma::uvec remaining = setdiff(indices, target);

  arma::mat Xcollapsible = X.cols(target);
  arma::vec Wcollapsible = w(target);

  arma::uvec weight_sort = arma::sort_index(Wcollapsible);

  for(arma::uword i = 0; i < weight_sort.n_rows / 5; i++) {
    if(i == weight_sort.n_rows / 5 - 1 || i * 5 + 4 >= weight_sort.n_rows) {
      arma::span cur_span(i * 5, weight_sort.n_rows - 1);
      arma::mat Xmat = Xcollapsible.cols(weight_sort(cur_span));
      arma::vec weight_spanned = Wcollapsible(weight_sort(cur_span));
      groups.push_back(VariantGroup(Xmat, Y, weight_spanned, group_threshold, site_penalty));
    } else {
	  arma::span cur_span(i * 5, i * 5 + 4);
	  arma::mat Xmat = Xcollapsible.cols(weight_sort(cur_span));
	  arma::vec weight_spanned = Wcollapsible(weight_sort(cur_span));
	  groups.push_back(VariantGroup(Xmat, Y, weight_spanned, group_threshold, site_penalty));
    }
  }
  // Don't group the remaining variants.
  if(remaining.n_rows > 0) {
	groups.push_back(VariantGroup(X.cols(remaining), Y, w.rows(remaining), 0, site_penalty));
  }
}

void VAAST::variant_bitmask(const arma::mat &X, const arma::vec &Y, const arma::vec &w) {
  // mask variants with score < 1
  arma::uvec tmpmask = arma::find(vaast_site_scores <= 0);
  std::vector<arma::uword> scenario;
  /*
   * Collect group variants.
   */
  arma::mat Xmat;

  for(auto it = groups.begin(); it != groups.end(); it++) {
    if(it == groups.begin()) {
      Xmat = (*it).X;
    } else {
	  Xmat = join_horiz(X, (*it).X);
    }
  }

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
  case_allele1 = Xmat.t() * Y;
  control_allele1 = Xmat.t() * (1 - Y);
  case_allele0 = 2. * n_case - case_allele1;
  control_allele0 = 2. * n_control - control_allele1;

  for (arma::uword i = 0; i < case_allele1.n_rows; i++) {
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
