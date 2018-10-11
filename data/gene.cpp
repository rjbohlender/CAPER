//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include <string>

#include "gene.hpp"
#include <boost/format.hpp>

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

arma::mat &Gene::get_matrix(const std::string &k) {
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

void Gene::clear(Covariates &cov, std::unordered_map<std::string, Result> &results) {
  generate_detail(cov, results);
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
	  for(arma::uword j = 3; j < splitter.size(); j++) {
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
	  genotypes_[transcripts_.back()] = arma::mat(nsamples_, nvariants_[transcripts_.back()]);
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
	  try{
		val = std::stod(splitter[j]);
	  } catch(std::exception &e) {
	    std::cerr << "Full buffer: " << ss.str() << std::endl;
	    std::cerr << "Failed to convert data to double: " << splitter[j] << std::endl;
	    std::cerr << "Line: " << line << std::endl;
	    std::cerr << "j: " << j << std::endl;
	    std::exit(-1);
	  }
	  // Handle missing data
	  if (val > 2 || val < 0)
		val = 0;
	  genotypes_[transcripts_.back()](j - 3, i - 1) = val;
	}
	i++;
  }
  // Switch to counting minor allele
  for (auto &v : genotypes_) {
	// For each variant
	for (i = 0; i < v.second.n_cols; i++) {
	  // Check allele frequency
	  if (arma::mean(v.second.col(i)) / 2 > 0.5) {
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
	if(weights.n_rows != genotypes_[k].n_cols) {
	  throw(std::logic_error("Weights do not match number of variants."));
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

void Gene::generate_detail(Covariates &cov, std::unordered_map<std::string, Result> &results) {
  std::stringstream detail;

  arma::vec Y = cov.get_original_phenotypes();

  std::map<std::string, std::vector<std::string>> pos_ts_map;
  std::map<std::string, double> pos_score_map;
  std::map<std::string, double> pos_freq_map;

  // Case/Control Ref and Alt counts
  std::map<std::string, double> pos_caseref_map;
  std::map<std::string, double> pos_casealt_map;
  std::map<std::string, double> pos_contref_map;
  std::map<std::string, double> pos_contalt_map;

  // Alt carriers
  std::map<std::string, arma::uvec> pos_caseidx_map;
  std::map<std::string, arma::uvec> pos_contidx_map;

  // Collect all positions across transcripts and associated scores
  for(const auto &ts : transcripts_) {
    arma::uword i = 0;

    arma::uvec cases = arma::find(Y == 1);
    arma::uvec controls = arma::find(Y == 0);

    arma::mat X = genotypes_[ts];
    arma::mat Xcase = X.rows(cases);
    arma::mat Xcont = X.rows(controls);
    arma::rowvec maf = arma::mean(X, 0) / 2.;
    // Ref/Alt Counts
    arma::rowvec case_alt = arma::sum(Xcase, 0);
	arma::rowvec case_ref = 2 * Xcase.n_rows - case_alt;
	arma::rowvec cont_alt = arma::sum(Xcont, 0);
	arma::rowvec cont_ref = 2 * Xcont.n_rows - cont_alt;

    for(const auto &pos : positions_[ts]) {
      // Get transcripts
      if(pos_ts_map.find(pos) == pos_ts_map.end()) {
        pos_ts_map[pos] = {ts};
      } else {
        pos_ts_map[pos].push_back(ts);
      }
      // Get scores
      if(pos_score_map.find(pos) == pos_score_map.end()) {
        if(variant_scores_[ts].empty())
          variant_scores_[ts].zeros(positions_[ts].size());
        pos_score_map[pos] = variant_scores_[ts](i);
        results[ts].testable = arma::find(variant_scores_[ts] > 0).eval().n_elem >= 4;
      }
      // Get frequency
      if(pos_freq_map.find(pos) == pos_freq_map.end()) {
        pos_freq_map[pos] = maf(i);
      }
      // Get counts
	  if(pos_caseref_map.find(pos) == pos_caseref_map.end()) {
		pos_caseref_map[pos] = case_ref(i);
	  }
	  if(pos_casealt_map.find(pos) == pos_casealt_map.end()) {
		pos_casealt_map[pos] = case_alt(i);
	  }
	  if(pos_contref_map.find(pos) == pos_contref_map.end()) {
		pos_contref_map[pos] = cont_ref(i);
	  }
	  if(pos_contalt_map.find(pos) == pos_contalt_map.end()) {
		pos_contalt_map[pos] = cont_alt(i);
	  }
	  // Get indices
	  arma::uvec carriers = arma::find(X.col(i) > 0);
	  if(pos_caseidx_map.find(pos) == pos_caseidx_map.end()) {
	    pos_caseidx_map[pos] = arma::intersect(cases, carriers);
	  }
	  if(pos_contidx_map.find(pos) == pos_contidx_map.end()) {
		pos_contidx_map[pos] = arma::intersect(controls, carriers);
	  }
      i++;
    }
  }

  for(const auto &v : pos_ts_map) {
    detail << gene_ << "\t";
    print_comma_sep(v.second, detail);
    detail << "\t";
    detail << boost::format("%1$s\t%2$.2f\t%3$d\t%4$d\t%5$d\t%6$d\t%7$d\t")
    	% v.first
    	% pos_score_map[v.first]
    	% pos_freq_map[v.first]
    	% pos_caseref_map[v.first]
    	% pos_casealt_map[v.first]
    	% pos_contref_map[v.first]
    	% pos_contalt_map[v.first];
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

void print_comma_sep(arma::uvec &x, std::ostream &os) {
  for(arma::uword i = 0; i < x.size(); i++) {
	if(i == x.size() - 1) {
	  os << x(i);
	} else {
	  os << x(i) << ",";
	}
  }
}

void print_comma_sep(std::vector<std::string> &x, std::ostream &os) {
  if(x.size() == 0) {
	os << "-";
	return;
  }
  for(arma::uword i = 0; i < x.size(); i++) {
	if(i == x.size() - 1) {
	  os << x.at(i);
	} else {
	  os << x.at(i) << ",";
	}
  }
}

void print_comma_sep(const std::vector<std::string> &x, std::ostream &os) {
  if(x.size() == 0) {
    os << "-";
    return;
  }
  for(arma::uword i = 0; i < x.size(); i++) {
	if(i == x.size() - 1) {
	  os << x.at(i);
	} else {
	  os << x.at(i) << ",";
	}
  }
}
void print_semicolon_sep(arma::uvec &x, std::ostream &os) {
  if(x.size() == 0) {
	os << "-";
	return;
  }
  for(arma::uword i = 0; i < x.size(); i++) {
	if(i == x.size() - 1) {
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

