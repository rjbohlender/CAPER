//
// Created by Bohlender,Ryan James on 10/2/18.
//

#define PRINT_MERGEINFO 0
#include "vaast.hpp"

#include <algorithm>

VAASTLogic::VAASTLogic(Gene &gene,
					   arma::vec &Y_,
					   const std::string &k,
					   bool score_only_minor,
					   bool score_only_alternative,
					   double site_penalty,
					   arma::uword group_threshold,
					   bool detail,
					   bool biallelic,
					   double soft_maf_filter,
					   bool legacy_)
	: som(score_only_minor), soa(score_only_alternative), detail(detail), biallelic(biallelic), k(k),
	  group_threshold(group_threshold), soft_maf_filter(soft_maf_filter),
	  site_penalty(site_penalty), legacy(legacy_),
	  X(gene.get_matrix(k)), Y(Y_) {

  // Verify weights are okay
  check_weights(gene);

  n_case = arma::accu(Y);
  n_control = arma::accu(1. - Y);

  if (group_threshold == 0 || X.n_cols == 1) {
	score = Score(X, Y, weights);
	if (detail) {
	  gene.set_scores(k, vaast_site_scores);
	  expanded_scores = vaast_site_scores;
	}
  } else {
	score = Score(X, Y, weights);
	if(legacy) {
	  legacy_grouping(X, Y, weights, gene.get_positions(k));
	} else {
	  variant_grouping(X, Y, weights, gene.get_positions(k));
	}
	score = Score();
  }

  if (detail) {
	if (group_threshold == 0 || X.n_cols == 1) {
	  gene.set_scores(k, vaast_site_scores);
	} else {
	  expanded_scores.reshape(X.n_cols, 1);
	  if (X.n_cols > 1) {
		arma::sp_mat Xnew(X);
		arma::vec Wnew(weights);
		std::vector<std::string> positions = gene.get_positions(k);
		for (arma::uword i = 0; i < X.n_cols; i++) {
		  if (i == 0) {
			positions.erase(positions.begin());
			Xnew.shed_col(i);
			Wnew.shed_row(i);
		  } else {
			positions[i - 1] = gene.get_positions(k)[i - 1];
			// TODO Find a way to avoid this copy, the old way gave a spontaneous seg fault
			Xnew = X;
			Xnew.shed_col(i);
			Wnew(i - 1) = weights(i - 1);
		  }

		  Score(Xnew, Y, Wnew); // Do normal individual variant scoring
		  if(legacy) {
			legacy_grouping(Xnew, Y, Wnew, positions);
		  } else {
			variant_grouping(Xnew, Y, Wnew, positions);
		  }
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

VAASTLogic::VAASTLogic(arma::sp_mat X_,
					   arma::vec &Y_,
					   arma::vec &weights,
					   std::vector<std::string> &positions_,
					   std::string k,
					   bool score_only_minor,
					   bool score_only_alternative,
					   bool biallelic,
					   arma::uword group_threshold,
					   double site_penalty,
					   double soft_maf_filter,
					   bool legacy_)
	: som(score_only_minor), soa(score_only_alternative), detail(true), biallelic(biallelic), k(std::move(k)),
	  group_threshold(group_threshold), soft_maf_filter(soft_maf_filter),
	  site_penalty(site_penalty), legacy(legacy_),
	  X(std::move(X_)), Y(Y_) {

  if (group_threshold == 0 || X.n_cols == 1) {
	score = Score(X, Y, weights);
	expanded_scores = vaast_site_scores;
  } else {
	score = Score(X, Y, weights);
	if(legacy) {
	  legacy_grouping(X, Y, weights, positions_);
	} else {
	  variant_grouping(X, Y, weights, positions_);
	}
	score = Score();
  }

  if (group_threshold > 0 && X.n_cols > 1) {
	expanded_scores.reshape(X.n_cols, 1);
	arma::sp_mat Xnew(X);
	arma::vec Wnew(weights);
	std::vector<std::string> positions = positions_;
	for (arma::uword i = 0; i < X.n_cols; i++) {
	  if (i == 0) {
		positions.erase(positions.begin());
		Xnew.shed_col(i);
		Wnew.shed_row(i);
	  } else {
		positions[i - 1] = positions_[i - 1];
		Xnew.col(i - 1) = X.col(i - 1);
		Wnew(i - 1) = weights(i - 1);
	  }

	  Score(Xnew, Y, Wnew); // Do normal individual variant scoring
	  if(legacy) {
		legacy_grouping(Xnew, Y, Wnew, positions);
	  } else {
		variant_grouping(Xnew, Y, Wnew, positions);
	  }
	  double new_score = Score();
	  expanded_scores(i) = (score - new_score >= 0) ? score - new_score : 0; // Calculate score
	}
  } else {
	expanded_scores = vaast_site_scores;
  }
}

double VAASTLogic::Score() {
  double val = 0;
  for (auto &s : vaast_site_scores) {
	if (s > site_penalty) {
	  s -= site_penalty;
	  val += s;
	}
  }
  return val;
}

double VAASTLogic::Score(const arma::sp_mat &X, const arma::vec &Y, const arma::vec &w) {
  // Get sample sizes
  n_case = static_cast<arma::uword>(arma::sum(Y));
  n_control = static_cast<arma::uword>(arma::sum(1 - Y));

  // Get carrier counts
  case_allele1 = X.t() * Y;
  control_allele1 = X.t() * (1 - Y);
  case_allele0 = 2. * n_case - case_allele1;
  control_allele0 = 2. * n_control - control_allele1;

  // Handle collapsed variant
  if (biallelic) {
	// Via email:
	// The idea is to add an additional term to the VAAST CLRT that captures the frequency of biallelic variants in cases and controls.
	// For starters, we can try collapsing all variants with a positive VAAST score into one biallelic term.
	arma::vec log_lh = LRT();
	vaast_site_scores = 2.0 * (log_lh + arma::log(w));
	variant_bitmask(X, Y, w);
	vaast_site_scores(mask).zeros();

	vaast_site_scores(arma::find(vaast_site_scores <= 2)).zeros();
	vaast_site_scores(arma::find(vaast_site_scores > 2)) -= site_penalty;

	// Get all variants with positive vaast score
	arma::mat Xmat(X);
	arma::uvec collapse = arma::find(vaast_site_scores > 0);
	arma::vec variant(n_case + n_control, arma::fill::zeros);

	variant = arma::sum(Xmat.cols(collapse), 1);
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

  arma::uvec unweighted = arma::find(w == 0);
  arma::vec log_lh = LRT();
  vaast_site_scores = 2.0 * (log_lh + arma::log(w));
  vaast_site_scores(arma::find(vaast_site_scores / (-2.) > 0)).zeros();
  vaast_site_scores(unweighted).zeros();
  variant_bitmask(X, Y, w); // Mask variants with CASM weight 0
  vaast_site_scores(mask).zeros();

  arma::uvec included = arma::find(vaast_site_scores > 2);
  vaast_site_scores(included) -= site_penalty;

  // vaast_site_scores(arma::find(vaast_site_scores <= 2)).zeros();
  double val = arma::accu(vaast_site_scores(included));
  // val += arma::accu(vaast_site_scores(arma::find(vaast_site_scores <= 2)));
  return (val >= 0) ? val : 0; // Mask all negative values
}

arma::vec VAASTLogic::LRT() {
  arma::vec alt_control_freq = control_allele1 / (control_allele0 + control_allele1);
  arma::vec alt_case_freq = case_allele1 / (case_allele0 + case_allele1);

  // Match -r behavior for VAAST
  alt_control_freq(arma::find(alt_control_freq > soft_maf_filter)).fill(soft_maf_filter);

  arma::vec
	  null_freq = (case_allele1 + control_allele1) / (case_allele0 + case_allele1 + control_allele0 + control_allele1);

  arma::vec alt_control = log_likelihood(alt_control_freq, control_allele0, control_allele1);
  arma::vec alt_case = log_likelihood(alt_case_freq, case_allele0, case_allele1);
  alt_control.replace(arma::datum::nan, 0);

  arma::vec alt_log_lh = alt_control + alt_case;

  arma::vec null_control = log_likelihood(null_freq, control_allele0, control_allele1);
  arma::vec null_case = log_likelihood(null_freq, case_allele0, case_allele1);
  arma::vec null_log_lh = null_case + null_control;

  arma::vec log_lh = alt_log_lh - null_log_lh;

  log_lh(arma::find(alt_case_freq == 0)).zeros();
  return log_lh;
}

arma::vec VAASTLogic::log_likelihood(arma::vec &freq, arma::vec &allele0, arma::vec &allele1) {
  // Prevent numerical issues
  // arma::vec clamped = arma::clamp(freq, 1e-9, 1.0 - 1e-9);

  return allele1 % arma::log(freq) + allele0 % arma::log(1.0 - freq);
}

void VAASTLogic::check_weights(Gene &gene) {
  if (gene.is_weighted(k)) {
	weights = gene.get_weights(k);
	return;
  }

  weights = arma::vec(X.n_cols, arma::fill::ones);
  gene.set_weights(k, weights);
}

//! \brief Calculate the set difference between two uvec
//! \param x uvec with elements to keep
//! \param y uvec with elements to exclude
//! \return copy of uvec excluding any element in y
arma::uvec VAASTLogic::setdiff(arma::uvec x, arma::uvec y) {
  for (size_t j = 0; j < y.n_elem; j++) {
	arma::uvec q1 = arma::find(x == y[j]);
	if (!q1.empty()) {
	  x.shed_row(q1(0));
	}
  }

  return x;
}

void VAASTLogic::variant_bitmask(const arma::sp_mat &X, const arma::vec &Y, const arma::vec &w) {
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

void VAASTLogic::variant_grouping(const arma::sp_mat &X,
								  const arma::vec &Y,
								  const arma::vec &w,
								  std::vector<std::string> &positions) {
  case_allele1 = X.t() * Y;
  control_allele1 = X.t() * (1 - Y);
  case_allele0 = 2 * n_case - case_allele1;
  control_allele0 = 2 * n_control - control_allele1;

  std::vector<Variant> variants;

  for (int i = 0; i < positions.size(); i++) {
    RJBUtil::Splitter<std::string> splitter(positions[i], "-");
	std::stringstream loc;
	loc << splitter[0] << "-" << splitter[1] << "-" << splitter[2];
	variants.emplace_back(Variant(case_allele1(i),
								  case_allele0(i),
								  control_allele1(i),
								  control_allele0(i),
								  w(i),
								  splitter.back(),
								  loc.str(),
								  i,
								  soft_maf_filter));
  }

  std::vector<std::vector<Variant>> sub_groups;
  std::vector<Variant> unmerged;
  std::vector<Variant> grouped;

#if PRINT_MERGEINFO
  if(!printed_mergeinfo) {
	std::stringstream ss;
	ss << "/Users/rjbohlender/CLionProjects/PA_Test/aas_no_intron/" << k << ".mergeinfo";
	merge_debug.open(ss.str());
  }
#endif

  std::stable_sort(variants.begin(), variants.end(), [&](auto &a, auto &b) {
	if (w(a.index) == w(b.index)) {
	  if (a.score == b.score) {
		return positions[a.index] > positions[b.index]; // Lexicographically sorted to match default VAAST behavior
	  } else {
		return a.score > b.score;
	  }
	} else {
	  return w(a.index) > w(b.index);
	}
  });

  int skipped = 0;
  std::vector<Variant> merging_vars;
  for (const auto &v : variants) {
	if (v.case_allele1 > group_threshold || v.score <= 0) {
	  unmerged.push_back(v);
	  skipped++;
	} else {
	  merging_vars.push_back(v);
	}
  }

  int sites = merging_vars.size() / 5;
  for (int i = 0; i < merging_vars.size(); i++) {
	if (sub_groups.size() < sites || sub_groups.size() == 0) {
	  if ((i % 5) == 0) {
		sub_groups.emplace_back(std::vector<Variant>());
	  }
	}
	// std::cerr << "gid: " << sub_groups.size() << " " << "var: " << positions[type.second[i].index] << " score: " << type.second[i].score << std::endl;
	sub_groups.back().emplace_back(std::move(merging_vars[i]));
  }

  int gid = 0;
  for (auto &group : sub_groups) {
	gid++;
	// std::cerr << "gid: " << gid << " group_size:" << group.size() << std::endl;
	std::stable_sort(group.begin(), group.end(), [](auto &a, auto &b) { return a.score < b.score; });
	Variant old("SNV", soft_maf_filter);
	Variant merged("SNV", soft_maf_filter);
	old.merge(group[group.size() - 1]);
	merged.merge(group[group.size() - 1]);
#if PRINT_MERGEINFO
	if(!printed_mergeinfo)
	  merge_debug << gid << " " << positions[group.back().index] << std::endl;
#endif
	for (int j = static_cast<int>(group.size()) - 2; j >= 0; j--) {
	  merged.merge(group[j]);
	  if (merged.score > old.score) {
#if PRINT_MERGEINFO
		if (!printed_mergeinfo)
		  merge_debug << gid << " " << positions[group[j].index] << std::endl;
#endif
		old = merged;
	  } else {
		merged = old;
		unmerged.push_back(group[j]);
	  }
	}
	grouped.push_back(merged);
  }

  vaast_site_scores.resize(unmerged.size() + grouped.size());
  vaast_site_scores.zeros();
  int i = 0;
  for (auto &v : unmerged) {
#if PRINT_MERGEINFO
	if(!printed_mergeinfo)
	  merge_debug << "unmerged " << positions[v.index] << std::endl;
#endif
	vaast_site_scores(i) = v.score;
	i++;
  }
  for (auto &v : grouped) {
	vaast_site_scores(i) = v.score;
	i++;
  }
}

void VAASTLogic::legacy_grouping(const arma::sp_mat &X,
								 const arma::vec &Y,
								 const arma::vec &w,
								 std::vector<std::string> &positions) {
  case_allele1 = X.t() * Y;
  control_allele1 = X.t() * (1 - Y);
  case_allele0 = 2 * n_case - case_allele1;
  control_allele0 = 2 * n_control - control_allele1;

  std::map<std::string, std::vector<Variant>> type_idx_map;
  std::vector<Variant> SNVs;
  std::vector<Variant> insertion;
  std::vector<Variant> deletion;
  std::vector<Variant> mnp;
  std::vector<Variant> spda;

  for (int i = 0; i < positions.size(); i++) {
	RJBUtil::Splitter<std::string> splitter(positions[i], "-");
	std::stringstream loc;
	loc << splitter[0] << "-" << splitter[1] << "-" << splitter[2];
	if (splitter.back() == "SNV") {
	  SNVs.emplace_back(Variant(case_allele1(i),
								case_allele0(i),
								control_allele1(i),
								control_allele0(i),
								w(i),
								splitter.back(),
								loc.str(),
								i,
								soft_maf_filter));
	} else if (splitter.back() == "insertion") {
	  insertion.emplace_back(Variant(case_allele1(i),
									 case_allele0(i),
									 control_allele1(i),
									 control_allele0(i),
									 w(i),
									 splitter.back(),
									 loc.str(),
									 i,
									 soft_maf_filter));
	} else if (splitter.back() == "deletion") {
	  deletion.emplace_back(Variant(case_allele1(i),
									case_allele0(i),
									control_allele1(i),
									control_allele0(i),
									w(i),
									splitter.back(),
									loc.str(),
									i,
									soft_maf_filter));
	} else if (splitter.back() == "SPDA") {
	  deletion.emplace_back(Variant(case_allele1(i),
									case_allele0(i),
									control_allele1(i),
									control_allele0(i),
									w(i),
									splitter.back(),
									loc.str(),
									i,
									soft_maf_filter));
	} else if (splitter.back() == "complex_substitution") {
	  mnp.emplace_back(Variant(case_allele1(i),
							   case_allele0(i),
							   control_allele1(i),
							   control_allele0(i),
							   w(i),
							   splitter.back(),
							   loc.str(),
							   i,
							   soft_maf_filter));
	} else {
	  throw (std::logic_error("Wrong variant type in VAAST legacy grouping."));
	}
  }
  type_idx_map["SNV"] = std::move(SNVs);
  type_idx_map["insertion"] = std::move(insertion);
  type_idx_map["deletion"] = std::move(deletion);
  type_idx_map["mnp"] = std::move(mnp);
  type_idx_map["spda"] = std::move(spda);

  std::vector<Variant> unmerged;
  std::vector<Variant> grouped;

  std::ofstream merge_debug;
#if PRINT_MERGEINFO
  if(!printed_mergeinfo) {
	std::stringstream ss;
	ss << "/Users/rjbohlender/CLionProjects/PA_Test/aas_no_intron/" << k << ".mergeinfo";
	merge_debug.open(ss.str());
  }
#endif

  for (auto &type : type_idx_map) {
	std::vector<std::vector<Variant>> sub_groups;

	std::stable_sort(type.second.begin(), type.second.end(), [&](auto &a, auto &b) {
	  if (w(a.index) == w(b.index)) {
		if (a.score == b.score) {
		  return positions[a.index] > positions[b.index]; // Lexicographically sorted to match default VAAST behavior
		} else {
		  return a.score > b.score;
		}
	  } else {
		return w(a.index) > w(b.index);
	  }
	});

	int skipped = 0;
	std::vector<Variant> merging_vars;
	for (const auto &v : type.second) {
	  if (v.case_allele1 > group_threshold || v.score <= 0) {
		unmerged.push_back(v);
		skipped++;
	  } else {
		merging_vars.push_back(v);
	  }
	}

	int sites = merging_vars.size() / 5;
	for (int i = 0; i < merging_vars.size(); i++) {
	  if (sub_groups.size() < sites || sub_groups.size() == 0) {
		if ((i % 5) == 0) {
		  sub_groups.emplace_back(std::vector<Variant>());
		}
	  }
	  // std::cerr << "gid: " << sub_groups.size() << " " << "var: " << positions[type.second[i].index] << " score: " << type.second[i].score << std::endl;
	  sub_groups.back().emplace_back(std::move(merging_vars[i]));
	}

	int gid = 0;
	for (auto &group : sub_groups) {
	  gid++;
	  // std::cerr << "gid: " << gid << " group_size:" << group.size() << std::endl;
	  std::stable_sort(group.begin(), group.end(), [](auto &a, auto &b) { return a.score < b.score; });
	  Variant old(type.first, soft_maf_filter);
	  Variant merged(type.first, soft_maf_filter);
	  old.merge(group[group.size() - 1]);
	  merged.merge(group[group.size() - 1]);
#if PRINT_MERGEINFO
	  if(!printed_mergeinfo)
		merge_debug << gid << " " << positions[group.back().index] << std::endl;
#endif
	  for (int j = static_cast<int>(group.size()) - 2; j >= 0; j--) {
		merged.merge(group[j]);
		if (merged.score > old.score) {
		  if (!printed_mergeinfo)
			merge_debug << gid << " " << positions[group[j].index] << std::endl;
		  old = merged;
		} else {
		  merged = old;
		  unmerged.push_back(group[j]);
		}
	  }
	  grouped.push_back(merged);
	}
  }

  vaast_site_scores.resize(unmerged.size() + grouped.size());
  vaast_site_scores.zeros();
  int i = 0;
  for (auto &v : unmerged) {
#if PRINT_MERGEINFO
	if(!printed_mergeinfo)
	  merge_debug << "unmerged " << positions[v.index] << std::endl;
#endif
	vaast_site_scores(i) = v.score;
	i++;
  }
  for (auto &v : grouped) {
	vaast_site_scores(i) = v.score;
	i++;
  }
#if PRINT_MERGEINFO
  if(!printed_mergeinfo) {
	merge_debug.close();
	printed_mergeinfo = true;
  }
#endif
}

Variant::Variant(std::string type, double soft_maf_filter) : type(std::move(type)), soft_maf_filter(soft_maf_filter) {

}

Variant::Variant(double case_allele1,
				 double case_allele0,
				 double control_allele1,
				 double control_allele0,
				 double weight,
				 std::string type,
				 std::string loc,
				 int index,
				 double soft_maf_filter)
	: case_allele1(case_allele1), case_allele0(case_allele0), control_allele1(control_allele1),
	  control_allele0(control_allele0), weight(weight), type(std::move(type)), loc(std::move(loc)), index(index),
	  soft_maf_filter(soft_maf_filter) {
  calc_score();
}

void Variant::merge(Variant &other) {
  case_allele1 += other.case_allele1;
  case_allele0 += other.case_allele0;
  control_allele1 += other.control_allele1;
  control_allele0 += other.control_allele0;
  weight *= other.weight;

  calc_score();
}

void Variant::calc_score() {
  if (case_allele1 == 0) { // Background unique variants aren't scored.
	score = 0;
	return;
  }

  if (case_allele0 * control_allele1 > control_allele0 * case_allele1) { // OR filter in VAAST
	score = 0;
	return;
  }

  double case_alt_freq = case_allele1 / (case_allele0 + case_allele1);
  double control_alt_freq = control_allele1 / (control_allele0 + control_allele1);

  control_alt_freq = control_alt_freq > soft_maf_filter ? soft_maf_filter : control_alt_freq;

  double
	  null_freq = (case_allele1 + control_allele1) / (case_allele1 + case_allele0 + control_allele1 + control_allele0);

  double p_alt_case = case_allele1 * std::log(case_alt_freq) + case_allele0 * std::log(1. - case_alt_freq);
  double
	  p_alt_control = control_allele1 * std::log(control_alt_freq) + control_allele0 * std::log(1. - control_alt_freq);
  double p_null_case = case_allele1 * std::log(null_freq) + case_allele0 * std::log(1. - null_freq);
  double p_null_control = control_allele1 * std::log(null_freq) + control_allele0 * std::log(1. - null_freq);

  if (std::isnan(p_alt_case))
	p_alt_case = 0;
  if (std::isnan(p_alt_control))
	p_alt_control = 0;
  if (std::isnan(p_null_case))
	p_null_case = 0;
  if (std::isnan(p_null_control))
	p_null_control = 0;

  double alt_log_lik = p_alt_case + p_alt_control;
  double null_log_lik = p_null_case + p_null_control;

  if (std::isnan(alt_log_lik)) {
	alt_log_lik = 0;
  }
  if (std::isnan(null_log_lik)) {
	null_log_lik = 0;
  }

  double lrt = alt_log_lik - null_log_lik;

  if (weight <= 0) {
	score = 0;
  } else {
	score = 2 * (lrt + std::log(weight));
  }
  if (score / (-2.) > 0) {
	score = 0;
  }
}

