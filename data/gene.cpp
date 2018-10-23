//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include <string>
#include <unordered_map>

#include <boost/format.hpp>

#include "gene.hpp"
#include "../statistics/vaast.hpp"
#include "../link/binomial.hpp"
#include "../statistics/bayesianglm.hpp"

Gene::Gene(std::stringstream &ss,
		   unsigned long nsamples,
		   std::map<std::string, arma::uword> &nvariants,
		   const Weight &weight)
	: nsamples_(nsamples),
	  nvariants_(nvariants) {
  parse(ss);
  if (!weight.empty()) {
	for (const auto &k : transcripts_) {
	  weights_[k].reshape(nvariants_[k], 1);
	  arma::uword i = 0;
	  for (const auto &v : positions_[k]) {
		weights_[k](i) = weight.get(v);
		i++;
	  }
	  weights_set_[k] = true;
	}
  } else {
	for (const auto &k : transcripts_) {
	  weights_set_[k] = false;
	}
  }
}

void Gene::print() {
  for (const auto &v : genotypes_)
	std::cout << v.second;
}

arma::sp_mat &Gene::get_matrix(const std::string &k) {
  return genotypes_[k];
}

std::string &Gene::get_gene() {
  return gene_;
}

std::vector<std::string> &Gene::get_transcripts() {
  return transcripts_;
}

arma::vec &Gene::get_weights(const std::string &k) {
  return weights_[k];
}

arma::uword Gene::get_nvariants(const std::string &k) {
  return nvariants_[k];
}

std::vector<std::string> &Gene::get_positions(const std::string &k) {
  return positions_[k];
}

void Gene::clear(Covariates &cov, std::unordered_map<std::string, Result> &results, TaskParams &tp) {
  generate_detail(cov, results, tp);

  // Set matrix size to 0x0 to free space.
  for (auto &v : genotypes_) {
	v.second.reset();
  }
}

void Gene::parse(std::stringstream &ss) {
  std::string line;
  arma::uword i = 0;

  while (std::getline(ss, line, '\n')) {
	if (i == 0) {
	  header_ = line;
	  RJBUtil::Splitter<std::string> splitter(line, "\t");
	  for (arma::uword j = 3; j < splitter.size(); j++) {
		samples_.push_back(splitter[j]);
	  }
	  i++;
	  continue;
	}
	RJBUtil::Splitter<std::string> splitter(line, "\t");

	if (gene_.empty()) {
	  gene_ = splitter[0];
	}
	auto found = std::find(std::begin(transcripts_), std::end(transcripts_), splitter[1]);
	if (found == std::end(transcripts_)) {
	  // Transcript not found -- add
	  transcripts_.push_back(splitter[1]);
	  // Start with matrix transposed
	  genotypes_[transcripts_.back()] = arma::sp_mat(nsamples_, nvariants_[transcripts_.back()]);
	  weights_[transcripts_.back()] = arma::vec(nvariants_[transcripts_.back()], arma::fill::zeros);
	  // Reset counter on new transcript
	  i = 1;
	}
	if (positions_.find(transcripts_.back()) == positions_.end()) {
	  positions_[transcripts_.back()] = std::vector<std::string>();
	}
	positions_[transcripts_.back()].push_back(splitter[2]);

	for (arma::uword j = 3; j < splitter.size(); j++) {
	  auto val = -1.;
	  try {
		val = std::stod(splitter[j]);
	  } catch (std::exception &e) {
		std::cerr << "Full buffer: " << ss.str() << std::endl;
		std::cerr << "Failed to convert data to double: " << splitter[j] << std::endl;
		std::cerr << "Line: " << line << std::endl;
		std::cerr << "j: " << j << std::endl;
		std::exit(-1);
	  }
	  // Handle missing data
	  if (val > 2 || val < 0)
		val = 0;
	  if(val != 0)
		genotypes_[transcripts_.back()](j - 3, i - 1) = val;
	}
	i++;
  }
  // Switch to counting minor allele
  for (auto &v : genotypes_) {
	// For each variant
	for (i = 0; i < v.second.n_cols; i++) {
	  // Check allele frequency
	  if (arma::mean(arma::vec(v.second.col(i))) / 2 > 0.5) {
		for (arma::uword j = 0; j < v.second.n_rows; j++) {
		  switch ((int) v.second(j, i)) {
		  case 0: {
			v.second(j, i) = 2;
			break;
		  }
		  case 1: {
			break;
		  }
		  case 2: {
			v.second(j, i) = 0;
			break;
		  }
		  default: break;
		  }
		}
	  }
	}
  }
}

void Gene::set_weights(const std::string &k, arma::vec &weights) {
  if (weights.n_rows != genotypes_[k].n_cols) {
	throw (std::logic_error("Weights do not match number of variants."));
  }

  weights_[k] = weights;

  weights_set_[k] = true;
}

bool Gene::is_weighted(const std::string &k) {
  return weights_set_[k];
}

arma::vec &Gene::get_scores(const std::string &k) {
  return variant_scores_[k];
}

void Gene::set_scores(const std::string &k, arma::vec &scores) {
  variant_scores_[k] = scores;
}

void Gene::generate_detail(Covariates &cov, std::unordered_map<std::string, Result> &results, TaskParams &tp) {
  std::stringstream detail;

  arma::vec Y = cov.get_original_phenotypes();

  std::unordered_map<std::string, std::vector<std::string>> pos_ts_map;
  std::unordered_map<std::string, double> pos_score_map;
  std::unordered_map<std::string, double> pos_freq_map;
  std::unordered_map<std::string, double> pos_odds_map;
  std::unordered_map<std::string, double> pos_odds_pval_map;

  // Case/Control Ref and Alt counts
  std::unordered_map<std::string, double> pos_caseref_map;
  std::unordered_map<std::string, double> pos_casealt_map;
  std::unordered_map<std::string, double> pos_contref_map;
  std::unordered_map<std::string, double> pos_contalt_map;

  // Alt carriers
  std::unordered_map<std::string, arma::uvec> pos_caseidx_map;
  std::unordered_map<std::string, arma::uvec> pos_contidx_map;

  // Collect all positions across transcripts and associated scores
  for (const auto &ts : transcripts_) {
	arma::uword i = 0;

	arma::uvec cases = arma::find(Y == 1);
	arma::uvec controls = arma::find(Y == 0);

	arma::sp_mat X(genotypes_[ts].t());
	arma::sp_mat Xcase(X.n_rows, cases.n_elem);
	arma::sp_mat Xcont(X.n_rows, controls.n_elem);

	arma::uword j = 0;
	for (const auto &k : cases) {
	  Xcase.col(j) = X.col(k);
	  j++;
	}
	j = 0;
	for (const auto &k : controls) {
	  Xcont.col(j) = X.col(k);
	  j++;
	}

	arma::vec maf = arma::vec(arma::mean(X, 1) / 2.);
	// Ref/Alt Counts
	arma::vec case_alt = X * Y;
	arma::vec case_ref = 2 * cases.n_elem - case_alt;
	arma::vec cont_alt = X * (1. - Y);
	arma::vec cont_ref = 2 * controls.n_elem - cont_alt;
	if (!tp.linear) {
	  // Get odds
	  Binomial link("logit");
	  arma::mat D = arma::join_vert(cov.get_covariate_matrix(), arma::mat(X).each_col() - arma::mean(arma::mat(X), 1));
	  BayesianGLM<Binomial> fit(D, Y, link);
	  for (const auto &pos : positions_[ts]) {
		// Get transcripts
		if (pos_ts_map.find(pos) == pos_ts_map.end()) {
		  pos_ts_map[pos] = {ts};
		} else {
		  pos_ts_map[pos].push_back(ts);
		}
		// Get scores
		if (pos_score_map.find(pos) == pos_score_map.end()) {
		  if (variant_scores_[ts].empty())
			variant_scores_[ts].zeros(positions_[ts].size());
		  pos_score_map[pos] = variant_scores_[ts](i);
		  pos_odds_map[pos] = fit.beta_(cov.get_covariate_matrix().n_rows + i);
		  pos_odds_pval_map[pos] = fit.pval_(cov.get_covariate_matrix().n_rows + i);
		}
		// Get frequency
		if (pos_freq_map.find(pos) == pos_freq_map.end()) {
		  pos_freq_map[pos] = maf(i);
		}
		// Get counts
		if (pos_caseref_map.find(pos) == pos_caseref_map.end()) {
		  pos_caseref_map[pos] = case_ref(i);
		}
		if (pos_casealt_map.find(pos) == pos_casealt_map.end()) {
		  pos_casealt_map[pos] = case_alt(i);
		}
		if (pos_contref_map.find(pos) == pos_contref_map.end()) {
		  pos_contref_map[pos] = cont_ref(i);
		}
		if (pos_contalt_map.find(pos) == pos_contalt_map.end()) {
		  pos_contalt_map[pos] = cont_alt(i);
		}
		// Get indices
		arma::uvec carriers = arma::find(arma::rowvec(X.row(i)) > 0);
		if (pos_caseidx_map.find(pos) == pos_caseidx_map.end()) {
		  pos_caseidx_map[pos] = arma::intersect(cases, carriers);
		}
		if (pos_contidx_map.find(pos) == pos_contidx_map.end()) {
		  pos_contidx_map[pos] = arma::intersect(controls, carriers);
		}
		i++;
	  }
	  results[ts].testable = testable(ts, cov, tp);
	} else {
	  for (const auto &pos : positions_[ts]) {
		// Get transcripts
		if (pos_ts_map.find(pos) == pos_ts_map.end()) {
		  pos_ts_map[pos] = {ts};
		} else {
		  pos_ts_map[pos].push_back(ts);
		}
		// Get scores
		if (pos_score_map.find(pos) == pos_score_map.end()) {
		  if (variant_scores_[ts].empty())
			variant_scores_[ts].zeros(positions_[ts].size());
		  pos_score_map[pos] = variant_scores_[ts](i);
		}
		// Get frequency
		if (pos_freq_map.find(pos) == pos_freq_map.end()) {
		  pos_freq_map[pos] = maf(i);
		}
		// Get counts
		if (pos_caseref_map.find(pos) == pos_caseref_map.end()) {
		  pos_caseref_map[pos] = case_ref(i);
		}
		if (pos_casealt_map.find(pos) == pos_casealt_map.end()) {
		  pos_casealt_map[pos] = case_alt(i);
		}
		if (pos_contref_map.find(pos) == pos_contref_map.end()) {
		  pos_contref_map[pos] = cont_ref(i);
		}
		if (pos_contalt_map.find(pos) == pos_contalt_map.end()) {
		  pos_contalt_map[pos] = cont_alt(i);
		}
		// Get indices
		arma::uvec carriers = arma::find(arma::rowvec(X.row(i)) > 0);
		if (pos_caseidx_map.find(pos) == pos_caseidx_map.end()) {
		  pos_caseidx_map[pos] = arma::intersect(cases, carriers);
		}
		if (pos_contidx_map.find(pos) == pos_contidx_map.end()) {
		  pos_contidx_map[pos] = arma::intersect(controls, carriers);
		}
		i++;
	  }
	}
  }

  for (const auto &v : pos_ts_map) {
	detail << gene_ << "\t";
	print_comma_sep(v.second, detail);
	detail << "\t";
	detail << boost::format("%1$s\t%2$.2f\t%3$.2f\t%4$.4f\t%5$d\t%6$d\t%7$d\t%8$d\t%9$d")
		% v.first
		% pos_score_map[v.first]
		% std::exp(pos_odds_map[v.first])
		% pos_odds_pval_map[v.first]
		% pos_freq_map[v.first]
		% pos_caseref_map[v.first]
		% pos_casealt_map[v.first]
		% pos_contref_map[v.first]
		% pos_contalt_map[v.first];
	detail << "\t";
	print_semicolon_sep(pos_caseidx_map[v.first], detail);
	detail << "\t";
	print_semicolon_sep(pos_contidx_map[v.first], detail);
	detail << std::endl;
  }
  detail_ = detail.str();
}

std::string Gene::get_detail() {
  return detail_;
}

std::vector<std::string> &Gene::get_samples() {
  return samples_;
}

auto Gene::testable(const std::string &k, Covariates &cov, TaskParams &tp) -> bool {
  arma::vec Yvec = cov.get_original_phenotypes();
  arma::sp_mat Xmat(genotypes_[k]);

  auto ncase = static_cast<arma::uword>(arma::accu(Yvec));

  // Most Extreme Phenotype Distribution
  arma::uvec mac_carriers = arma::vec(arma::sum(Xmat, 1)) > 0;
  arma::uvec maj_carriers = 1 - mac_carriers;

  arma::uword nmac = std::min(ncase, arma::accu(mac_carriers));
  arma::uword nmaj = (ncase - nmac >= 0) ? ncase - nmac : 0;

  mac_carriers = arma::find(arma::vec(arma::sum(Xmat, 1)) > 0);
  maj_carriers = arma::find(arma::vec(arma::sum(Xmat, 1)) == 0);

  arma::vec extreme_phen(Yvec.n_elem, arma::fill::zeros);
  extreme_phen(mac_carriers(arma::span(0, nmac - 1))).ones();
  if(nmaj > 0) {
    // There are more cases to distribute than minor allele carriers
    extreme_phen(maj_carriers(arma::span(0, nmaj - 1))).ones();
  }
  assert(arma::accu(extreme_phen) == ncase);

  VAAST vaast(genotypes_[k], cov, weights_[k], positions_[k], k, tp.score_only_minor, tp.score_only_alternative, 2., tp.group_size);

  return arma::accu(vaast.expanded_scores > 0) >= 4;
}

void print_comma_sep(arma::uvec &x, std::ostream &os) {
  for (arma::uword i = 0; i < x.size(); i++) {
	if (i == x.size() - 1) {
	  os << x(i);
	} else {
	  os << x(i) << ",";
	}
  }
}

void print_comma_sep(std::vector<std::string> &x, std::ostream &os) {
  if (x.empty()) {
	os << "-";
	return;
  }
  for (arma::uword i = 0; i < x.size(); i++) {
	if (i == x.size() - 1) {
	  os << x.at(i);
	} else {
	  os << x.at(i) << ",";
	}
  }
}

void print_comma_sep(const std::vector<std::string> &x, std::ostream &os) {
  if (x.empty()) {
	os << "-";
	return;
  }
  for (arma::uword i = 0; i < x.size(); i++) {
	if (i == x.size() - 1) {
	  os << x.at(i);
	} else {
	  os << x.at(i) << ",";
	}
  }
}
void print_semicolon_sep(arma::uvec &x, std::ostream &os) {
  if (x.empty()) {
	os << "-";
	return;
  }
  for (arma::uword i = 0; i < x.size(); i++) {
	if (i == x.size() - 1) {
	  os << x(i);
	} else {
	  os << x(i) << ";";
	}
  }
}

// Adapted from
// https://stackoverflow.com/questions/29724083/trying-to-write-a-setdiff-function-using-rcpparmadillo-gives-compilation-error
arma::uvec setdiff(arma::uvec x, arma::uvec y) {
  for (size_t j = 0; j < y.n_elem; j++) {
	arma::uvec q1 = arma::find(x == y[j]);
	if (!q1.empty()) {
	  x.shed_row(q1(0));
	}
  }

  return x;
}

