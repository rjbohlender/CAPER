//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include <string>
#include <unordered_map>

#include <boost/format.hpp>

#include "gene.hpp"
#include "../statistics/vaast.hpp"
#include "../link/binomial.hpp"
#include "../link/gaussian.hpp"
#include "../statistics/bayesianglm.hpp"
#include "../statistics/glm.hpp"
#include "../utility/indices.hpp"

Gene::Gene(std::stringstream &ss,
		   std::shared_ptr<Covariates> cov,
		   unsigned long nsamples,
		   std::map<std::string, arma::uword> &nvariants,
		   const Weight &weight,
		   TaskParams tp)
    : nsamples_(nsamples),
      nvariants_(nvariants),
      testable_(false),
      skippable_(false),
      tp_(std::move(tp)) {
  parse(ss, cov);
  if (!weight.empty()) {
    for (const auto &k : transcripts_) {
      weights_[k].reshape(nvariants_[k], 1);
      arma::uword i = 0;
      for (const auto &v : positions_[k]) {
        weights_[k](i) = weight.get(v);
        i++;
      }
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

void Gene::set_matrix(const std::string &k, arma::sp_mat &data) {
  genotypes_[k] = std::move(data);
}

void Gene::set_matrix(const std::string &k, arma::sp_mat &&data) {
  genotypes_[k] = std::move(data);
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

void Gene::clear(Covariates &cov, std::unordered_map<std::string, Result> &results, const TaskParams &tp) {
  if (!tp.no_detail)
    generate_detail(cov, results, tp);
  if (tp.method == "VAAST") {
    generate_vaast(cov);
  }

  // Set matrix size to 0x0 to free space.
  for (auto &v : genotypes_) {
    v.second.reset();
  }
}

void Gene::parse(std::stringstream &ss, std::shared_ptr<Covariates> cov) {
  std::string line;
  arma::uword i = 0;
  std::vector<bool> rarest;

  while (std::getline(ss, line, '\n')) {
    if (i == 0) {
	  RJBUtil::Splitter<std::string> splitter(line, "\t");
      for (arma::uword j = static_cast<int>(Indices::first); j < splitter.size(); j++) {
        samples_.push_back(splitter[j]);
      }
      i++;
      continue;
    }
    RJBUtil::Splitter<std::string> splitter(line, "\t");

    if (gene_name.empty()) {
	  gene_name = splitter[static_cast<int>(Indices::gene)];
    }
    auto iter = std::find(std::begin(transcripts_), std::end(transcripts_), splitter[static_cast<int>(Indices::transcript)]);
    std::string transcript = splitter[static_cast<int>(Indices::transcript)];
    if (iter == std::end(transcripts_)) {
      // Transcript not found -- add
      transcripts_.push_back(transcript);
      // Start with matrix transposed
      genotypes_.emplace(std::make_pair(transcript, arma::sp_mat(nsamples_, nvariants_[transcript])));
      // Separately record the controls
      missing_variant_carriers_.emplace(std::make_pair(transcript, arma::sp_mat(nsamples_, nvariants_[transcript])));
      // Ensure positions are allocated for current transcript
	  positions_[transcript] = std::vector<std::string>();
	  // Ensure weights are the right length
	  weights_[transcript].reshape(nvariants_[transcript], 1);
	  // Add transcript to the type map
	  type_.emplace(std::make_pair(transcript, std::vector<std::string>()));
	  // Add transcript to the reference allele map
	  reference_.emplace(std::make_pair(transcript, std::vector<std::string>()));
	  // Add transcript to the alternate allele map
	  alternate_.emplace(std::make_pair(transcript, std::vector<std::string>()));
	  // Add transcript to the variant function map
	  function_.emplace(std::make_pair(transcript, std::vector<std::string>()));
	  // Add transcript to the variant annotation map
	  annotation_.emplace(std::make_pair(transcript, std::vector<std::string>()));
	  // Reset counter on new transcript
      i = 1;
    }
    std::string var_id = form_variant_id(splitter);
	// Handle multiple copies of a variant e.g. triallelic. We only take the rarest.
    if(std::find(positions_[transcript].begin(), positions_[transcript].end(), var_id) != positions_[transcript].end()) {
    } else {
	  positions_[transcript].push_back(var_id);
	}
    weights_[transcript](i - 1) = std::stod(splitter[static_cast<int>(Indices::weight)]);

	// Record other info in line
	function_[transcript].push_back(splitter[static_cast<int>(Indices::function)]);
	annotation_[transcript].push_back(splitter[static_cast<int>(Indices::annotation)]);
	reference_[transcript].push_back(splitter[static_cast<int>(Indices::ref)]);
	alternate_[transcript].push_back(splitter[static_cast<int>(Indices::alt)]);
	type_[transcript].push_back(splitter[static_cast<int>(Indices::type)]);

	auto first_idx = static_cast<int>(Indices::first);
    for (auto j = first_idx; j < splitter.size(); j++) {
      double val;
      try {
        val = std::stod(splitter[j]);
      }
      catch (std::exception &e) {
        std::cerr << "Full buffer: " << ss.str() << std::endl;
        std::cerr << "Failed to convert data to double: " << splitter[j] << std::endl;
        std::cerr << "Line: " << line << std::endl;
        std::cerr << "j: " << j << std::endl;
        std::exit(-1);
      }
      // Handle missing data
      if (val > 2 || val < 0) {
        if (tp_.impute_to_mean) {
		  val = -9;
		} else {
		  val = 0;
		  missing_variant_carriers_[transcript](j - first_idx, i - 1) = 1;
		}
      }
      if (val != 0) {
        try {
          genotypes_[transcript](j - first_idx, i - 1) = val;
        } catch (std::exception &e) {
          std::cerr << "Failed to set value for: row = " << j - first_idx << " col = " << i - 1
                    << std::endl;
          std::cerr << "Gene: " << splitter[static_cast<int>(Indices::gene)]
          			<< " Transcript: " << splitter[static_cast<int>(Indices::transcript)]
                    << " Location: " << form_variant_id(splitter) << std::endl;
          throw (e);
        }
      }
    }
    i++;
  }
  if (tp_.impute_to_mean) {
	// Impute missing
	for (auto &v : genotypes_) {
	  arma::mat X(v.second);
	  arma::vec Y(cov->get_original_phenotypes());
	  arma::uvec cases = arma::find(Y == 1);
	  arma::uvec controls = arma::find(Y == 0);
	  for(arma::uword i = 0; i < X.n_cols; i++) {
		arma::vec Xc = X.col(i);
		arma::uvec nonmissing = arma::find(Xc != -9);
		arma::uvec missing = arma::find(Xc == -9);

		arma::uvec keep = arma::intersect(cases, nonmissing);
		if (keep.n_elem > 0) {
		  double maf = arma::sum(Xc(keep)) / (2. * keep.n_elem);
		  Xc(arma::intersect(cases, missing)).replace(-9, maf);
		}
		keep = arma::intersect(controls, nonmissing);
		if (keep.n_elem > 0) {
		  double maf = arma::sum(Xc(keep)) / (2. * keep.n_elem);
		  Xc(arma::intersect(controls, missing)).replace(-9, maf);
		}
		X.col(i) = Xc;
	  }
	  v.second = arma::sp_mat(X);
	}
  }
  // Switch to counting minor allele
  for (auto &v : genotypes_) {
    arma::rowvec maf = arma::rowvec(arma::mean(v.second) / 2.);
    // For each variant
    for (arma::uword k : arma::find(maf > 0.5).eval()) {
      // Swap assignment
      v.second.col(k).replace(0, 3); // Placeholder to unique value
      v.second.col(k).replace(2, 0);
      v.second.col(k).replace(3, 2);
    }
  }
  // Mark gene skippable to start
  skippable_ = true;
  // Drop variants with minor allele count > tp_.mac.
  for (const auto &ts : transcripts_) {
    arma::rowvec sums = arma::rowvec(arma::sum(genotypes_[ts], 0));

    arma::rowvec maf = arma::rowvec(arma::mean(genotypes_[ts]) / 2.);
    for (arma::sword k = sums.n_elem - 1; k >= 0; --k) {
      bool bmac = sums[k] > tp_.mac;
      bool bmaf = maf[k] > tp_.maf;
      if (bmac || bmaf) {
        if (tp_.verbose && bmac) {
          std::cerr << "Removing: " << gene_name << " " << ts << " " << positions_[ts][k] << " | count: "
					<< sums[k] << " due to MAC filter." << std::endl;
        } else if (tp_.verbose && bmaf) {
          std::cerr << "Removing: " << gene_name << " " << ts << " " << positions_[ts][k] << " | frequency: "
					<< maf[k] << " due to MAF filter." << std::endl;
        }
        sums.shed_col(k);
        genotypes_[ts].shed_col(k);
        weights_[ts].shed_row(k);
        missing_variant_carriers_[ts].shed_col(k);
        positions_[ts].erase(positions_[ts].begin() + k);
        function_[ts].erase(function_[ts].begin() + k);
		reference_[ts].erase(reference_[ts].begin() + k);
		alternate_[ts].erase(alternate_[ts].begin() + k);
		annotation_[ts].erase(annotation_[ts].begin() + k);
		type_[ts].erase(type_[ts].begin() + k);
        nvariants_[ts]--;
      }
    }
    // Check if any polymorphic. Mark transcripts skippable if all fixed.
    if (arma::accu(sums) < tp_.min_minor_allele_count || nvariants_[ts] < tp_.min_variant_count) {
      std::cerr << "gene: " << gene_name << " marked skippable." << std::endl;
      polymorphic_[ts] = false;
    } else {
      polymorphic_[ts] = true;
      skippable_ = false;
    }
  }
}

void Gene::set_weights(const std::string &k, arma::vec &weights) {
  if (weights.n_rows != genotypes_[k].n_cols) {
    throw (std::logic_error("Weights do not match number of variants."));
  }

  weights_[k] = weights;
}

arma::vec &Gene::get_scores(const std::string &k) {
  return variant_scores_[k];
}

void Gene::set_scores(const std::string &k, arma::vec &scores) {
  variant_scores_[k] = scores;
}

void Gene::generate_detail(Covariates &cov, std::unordered_map<std::string, Result> &results, const TaskParams &tp) {
  std::stringstream detail;

  arma::vec Y = cov.get_original_phenotypes();

  std::unordered_map<std::string, std::vector<std::string>> pos_ts_map;
  std::unordered_map<std::string, double> pos_score_map;
  std::unordered_map<std::string, double> pos_weight_map;
  std::unordered_map<std::string, double> pos_freq_map;
  std::unordered_map<std::string, double> pos_odds_map;
  std::unordered_map<std::string, double> pos_serr_map;
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

    arma::sp_mat X(genotypes_[ts]);

    arma::rowvec maf = arma::rowvec(arma::mean(X) / 2.);
    // Ref/Alt Counts
    arma::vec case_alt = X.t() * Y;
    arma::vec case_ref = 2 * cases.n_elem - case_alt;
    arma::vec cont_alt = X.t() * (1. - Y);
    arma::vec cont_ref = 2 * controls.n_elem - cont_alt;
    if (!tp.linear) {
      // Get odds
      Binomial link("logit");
      arma::mat D = arma::join_horiz(cov.get_covariate_matrix(),
                                     arma::mat(X).each_row() - 2 * maf);
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
          if (!fit.beta_.has_nan()) {
            pos_score_map[pos] = variant_scores_[ts](i);
            pos_weight_map[pos] = weights_[ts](i);
            pos_odds_map[pos] = fit.beta_(cov.get_covariate_matrix().n_cols + i);
            pos_serr_map[pos] = fit.s_err_(cov.get_covariate_matrix().n_cols + i);
            pos_odds_pval_map[pos] = fit.pval_(cov.get_covariate_matrix().n_cols + i);
          } else {
            // Regression failed -- too many features
            pos_score_map[pos] = variant_scores_[ts](i);
            pos_weight_map[pos] = weights_[ts](i);
            pos_odds_map[pos] = 1.;
            pos_serr_map[pos] = 1.;
            pos_odds_pval_map[pos] = 1.;
          }
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
        arma::uvec carriers = arma::find(arma::rowvec(X.col(i).t()) > 0);
        if (pos_caseidx_map.find(pos) == pos_caseidx_map.end()) {
          pos_caseidx_map[pos] = arma::intersect(cases, carriers);
        }
        if (pos_contidx_map.find(pos) == pos_contidx_map.end()) {
          pos_contidx_map[pos] = arma::intersect(controls, carriers);
        }
        i++;
      }
      if (tp.testable)
        results[ts].testable = testable(ts, cov, tp);
      if (!testable_)
        testable_ = results[ts].testable;
    } else {
      // Get odds via Moser & Coombs (2004) -- Rather than dichotomizing, we fit normally and recover OR
      Gaussian link("identity");
      arma::mat D = arma::join_vert(cov.get_covariate_matrix(),
                                    arma::mat(X).each_row() - maf);
      GLM<Gaussian> fit(D, Y, link, nullptr, TaskParams());

      // Transform values
      arma::uword n = cov.get_covariate_matrix().n_rows;
      double lambda = arma::datum::pi / std::sqrt(3);
      arma::vec var = link.variance(fit.mu_);
      arma::mat fisher_info = arma::inv(D * arma::diagmat(var) * D.t());
      arma::mat A = arma::eye(Y.n_elem, Y.n_elem) - D.t() * arma::inv(D * D.t()) * D;
      double sd = arma::as_scalar(arma::sqrt(Y.t() * A * Y / (Y.n_elem - D.n_rows)));

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
          pos_odds_map[pos] = lambda * fit.beta_(n + i - 1) / sd;
          pos_serr_map[pos] = std::sqrt(fisher_info.diag()(n + i));
          pos_odds_pval_map[pos] = arma::datum::nan;
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
  if (tp.linear) {
    // Output for quantitative traits
    for (const auto &v : pos_ts_map) {
      detail << gene_name << "\t";
      print_comma_sep(v.second, detail);
      detail << "\t";
      detail << boost::format("%1$s\t%2$.2f\t%3$.2f\t%4$d")
          % v.first
          % pos_score_map[v.first]
          % pos_weight_map[v.first]
          % pos_freq_map[v.first];
      detail << std::endl;
    }
  } else {
    // Output for binary traits
    for (const auto &v : pos_ts_map) {
      detail << gene_name << "\t";
      print_comma_sep(v.second, detail);
      detail << "\t";
      detail << boost::format("%1$s\t%2$.2f\t%3$.2f\t%4$d\t%5$d\t%6$d\t%7$d\t%8$d")
          % v.first
          % pos_score_map[v.first]
          % pos_weight_map[v.first]
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
  }
  detail_ = detail.str();
}

auto Gene::get_detail() const -> std::string {
  return detail_;
}

auto Gene::get_vaast() const -> std::map<std::string, std::string> {
  return vaast_;
}

std::vector<std::string> &Gene::get_samples() {
  return samples_;
}

auto Gene::testable(const std::string &k, Covariates &cov, const TaskParams &tp) -> bool {
  arma::vec Yvec = cov.get_original_phenotypes();
  arma::sp_mat Xmat(genotypes_[k]);

  auto ncase = static_cast<arma::sword>(arma::accu(Yvec));

  // Most Extreme Phenotype Distribution
  arma::uvec mac_carriers = arma::vec(arma::sum(Xmat, 1)) > 0;
  arma::uvec maj_carriers = 1 - mac_carriers;

  arma::sword nmac = std::min(ncase, static_cast<arma::sword>(arma::accu(mac_carriers)));
  arma::sword nmaj = std::max(ncase - nmac, 0ll);

  mac_carriers = arma::find(mac_carriers > 0);
  maj_carriers = arma::find(maj_carriers > 0);

  arma::vec extreme_phen(Yvec.n_elem, arma::fill::zeros);
  extreme_phen(mac_carriers(arma::span(0, nmac - 1))).ones();
  if (nmaj > 0) {
    // There are more cases to distribute than minor allele carriers
    extreme_phen(maj_carriers(arma::span(0, nmaj - 1))).ones();
  }
  assert(arma::accu(extreme_phen) == ncase);

  VAASTLogic vaast(genotypes_[k], extreme_phen, weights_[k], positions_[k], k, false, tp.group_size, 2., 0, false);

  return arma::accu(vaast.expanded_scores > 0) >= 4;
}

auto Gene::is_testable() const -> bool {
  return testable_;
}

auto Gene::is_skippable() const -> bool {
  return skippable_;
}

auto Gene::is_polymorphic(const std::string &k) -> bool {
  return polymorphic_[k];
}

auto Gene::generate_vaast(Covariates &cov) -> void {
  // One entry per transcript
  for (const auto &ts : transcripts_) {
    std::stringstream vaast_ss;
    vaast_ss << ">" << ts << "\t" << gene_name << std::endl;

    arma::mat X(genotypes_[ts]);
    arma::vec Y(cov.get_original_phenotypes());

    arma::uvec cases = arma::find(Y == 1);
    arma::uvec controls = arma::find(Y == 0);
    // Variant Format:
    // Variants are nucleotide then amino acid
    // R/U: variant_score position@chr nt|AA B/N|id1,id2...|ref:alt|ref:alt

    // Skip exon positions for now

    // Add target variants
    for (int i = 0; i < positions_[ts].size(); i++) {
      double case_alt = arma::as_scalar(X.col(i).t() * Y);
      double case_ref = 2 * cases.n_elem - case_alt;
      double cont_alt = arma::as_scalar(X.col(i).t() * (1. - Y));
      double cont_ref = 2 * controls.n_elem - cont_alt;

      if (case_alt > 0 && cont_alt == 0) {
        vaast_ss << "TU:\t";
      } else if (case_alt > 0 && cont_alt > 0) {
        vaast_ss << "TR:\t";
      } else {
        continue;
      }
      vaast_ss << boost::format("%1$.2f") % variant_scores_[ts][i] << "\t";

      RJBUtil::Splitter<std::string> pos_splitter(positions_[ts][i], "-");
      vaast_ss << pos_splitter[1] << "@" << pos_splitter[0] << "\t";

      // Placeholder for Reference alleles
      vaast_ss << boost::format("%1$s\t") % reference_[ts][i];
      vaast_ss << boost::format("%1$s\t") % annotation_[ts][i];

	  arma::uvec het_carriers = arma::find(X.col(i) % Y == 1);
	  arma::uvec hom_carriers = arma::find(X.col(i) % Y == 2);

	  if(het_carriers.n_elem > 0) {
		for (arma::uword j = 0; j < het_carriers.n_elem; j++) {
		  if (j < het_carriers.n_elem - 1) {
			vaast_ss << het_carriers(j) << ",";
		  } else {
			vaast_ss << het_carriers(j) << "|";
		  }
		}
		// Variant
		vaast_ss << boost::format("%1$s:%2$s") % reference_[ts][i] % alternate_[ts][i] << std::endl;
	  }
	  if(hom_carriers.n_elem > 0) {
		vaast_ss << "\t";
		for (arma::uword j = 0; j < hom_carriers.n_elem; j++) {
		  if (j < hom_carriers.n_elem - 1) {
			vaast_ss << hom_carriers(j) << ",";
		  } else {
			vaast_ss << hom_carriers(j) << "|";
		  }
		}
		// Variant
		vaast_ss << boost::format("%1$s:%2$s") % alternate_[ts][i] % alternate_[ts][i] << std::endl;
	  }
	}
    // Add background variants
    for (int i = 0; i < positions_[ts].size(); i++) {
      double case_alt = arma::as_scalar(X.col(i).t() * Y);
      double case_ref = 2 * cases.n_elem - case_alt;
      double cont_alt = arma::as_scalar(X.col(i).t() * (1. - Y));
      double cont_ref = 2 * controls.n_elem - cont_alt;

      if (case_alt == 0 && cont_alt > 0) {
        vaast_ss << "BU:\t";
      } else if (case_alt > 0 && cont_alt > 0) {
        vaast_ss << "BR:\t";
      } else {
        continue;
      }
      vaast_ss << boost::format("%1$.2f") % variant_scores_[ts][i] << "\t";

      RJBUtil::Splitter<std::string> pos_splitter(positions_[ts][i], "-");
      vaast_ss << pos_splitter[1] << "@" << pos_splitter[0] << "\t";

      // Placeholder for Reference alleles
	  vaast_ss << boost::format("%1$s\t") % reference_[ts][i];
	  vaast_ss << boost::format("%1$s\t") % annotation_[ts][i];

      arma::uvec het_carriers = arma::find(X.col(i) % (1. - Y) == 1);
	  arma::uvec hom_carriers = arma::find(X.col(i) % (1. - Y) == 2);

	  if(het_carriers.n_elem > 0) {
		for (arma::uword j = 0; j < het_carriers.n_elem; j++) {
		  if (j < het_carriers.n_elem - 1) {
			vaast_ss << het_carriers(j) << ",";
		  } else {
			vaast_ss << het_carriers(j) << "|";
		  }
		}
		// Variant
		vaast_ss << boost::format("%1$s:%2$s") % reference_[ts][i] % alternate_[ts][i] << std::endl;
	  }
	  if(hom_carriers.n_elem > 0) {
	    vaast_ss << "\t";
		for (arma::uword j = 0; j < hom_carriers.n_elem; j++) {
		  if (j < hom_carriers.n_elem - 1) {
			vaast_ss << hom_carriers(j) << ",";
		  } else {
			vaast_ss << hom_carriers(j) << "|";
		  }
		}
		// Variant
		vaast_ss << boost::format("%1$s:%2$s") % alternate_[ts][i] % alternate_[ts][i] << std::endl;
	  }
	}
	vaast_[ts] = vaast_ss.str();
  }
}

arma::sp_mat &Gene::get_missing(const std::string &k) {
  return missing_variant_carriers_[k];
}

std::string Gene::form_variant_id(RJBUtil::Splitter<std::string> &splitter) {
  std::stringstream ss;
  ss << splitter[static_cast<int>(Indices::chrom)] << "-"
  	 << splitter[static_cast<int>(Indices::start)] << "-"
  	 << splitter[static_cast<int>(Indices::end)] << "-"
  	 << splitter[static_cast<int>(Indices::type)];
  return ss.str();
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

