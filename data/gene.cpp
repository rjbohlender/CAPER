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

Gene::Gene(std::stringstream &ss, unsigned long nsamples, std::map<std::string, arma::uword> &nvariants,
           const Weight &weight, TaskParams tp, arma::vec &phenotypes)
    : nsamples_(nsamples),
      nvariants_(nvariants),
      tp_(std::move(tp)),
      testable_(false),
      skippable_(false) {
  parse(ss, phenotypes);
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

void Gene::set_matrix(const std::string &k, arma::sp_mat &&data) {
  genotypes_[k] = data;
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
  if (!tp.nodetail)
    generate_detail(cov, results, tp);
  if (tp.method == "VAAST") {
    generate_vaast(cov);
  }

  // Set matrix size to 0x0 to free space.
  for (auto &v : genotypes_) {
    v.second.reset();
  }
}

void Gene::parse(std::stringstream &ss, arma::vec &phenotypes) {
  std::string line;
  arma::uword i = 0;
  arma::uvec controls = arma::find(phenotypes == 0);
  arma::uword ncontrol = controls.n_elem;
  std::map<std::string, arma::sp_mat> control_gt;
  if (tp_.quantitative) { // Don't subset if not dichotomous
    controls = arma::conv_to<arma::uvec>::from(arma::regspace(0, 1, phenotypes.n_elem - 1));
    ncontrol = controls.n_elem;
  }

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
    std::string transcript = splitter[1];
    if (found == std::end(transcripts_)) {
      // Transcript not found -- add
      transcripts_.push_back(splitter[1]);
      // Start with matrix transposed
      genotypes_.emplace(std::make_pair(transcript, arma::sp_mat(nsamples_, nvariants_[transcript])));
      // Separately record the controls
      control_gt.emplace(std::make_pair(transcript, arma::sp_mat(ncontrol, nvariants_[transcript])));
      missing_variant_carriers_.emplace(std::make_pair(transcript, arma::sp_mat(nsamples_, nvariants_[transcript])));
      // Reset counter on new transcript
      i = 1;
    }
    if (positions_.find(transcript) == positions_.end()) {
      positions_[transcript] = std::vector<std::string>();
    }
    positions_[transcript].push_back(splitter[2]);

    arma::uword ncases = 0;
    for (arma::uword j = 3; j < splitter.size(); j++) {
      auto val = -1.;
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
#ifdef IMPUTE
        val = -9;
#else
        val = 0;
        missing_variant_carriers_[transcript](j - 3, i - 1) = 1;
#endif
      }
      if (val != 0) {
        try {
          genotypes_[transcript](j - 3, i - 1) = val;
          if(phenotypes(j - 3) == 0) { // Record control values
            control_gt[transcript](j - ncases - 3, i - 1) = val;
          } else {
            ncases++;
          }
        }
        catch (std::exception &e) {
          std::cerr << "Failed to set value for: row = " << j - 3 << " col = " << i - 1
                    << std::endl;
          std::cerr << "Gene: " << splitter[0] << " Transcript: " << splitter[1]
                    << " Location: " << splitter[2] << std::endl;
          throw (e);
        }
      }
    }
    i++;
  }
#ifdef IMPUTE
  // Impute missing
  for (auto &v : genotypes_) {
    arma::mat X(v.second);
    for(arma::uword i = 0; i < X.n_cols; i++) {
      arma::uvec nonmissing = arma::find(X.col(i) != -9);
      arma::vec Xc = X.col(i);

      double maf = arma::sum(Xc(nonmissing)) / (2. * nonmissing.n_elem);

      X.col(i).replace(-9, maf);
    }
    v.second = arma::sp_mat(X);
  }
#endif
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

    arma::sp_mat &gt_control = control_gt[ts];

    arma::rowvec maf = arma::rowvec(arma::mean(gt_control) / 2.);
    for (arma::sword i = sums.n_elem - 1; i >= 0; --i) {
      bool bmac = sums[i] > tp_.mac;
      bool bmaf = maf[i] > tp_.maf;
      if (bmac || bmaf) {
        if (tp_.verbose && bmac) {
          std::cerr << "Removing: " << gene_ << " " << ts << " " << positions_[ts][i] << " | count: "
                    << sums[i] << " due to MAC filter." << std::endl;
        } else if (tp_.verbose && bmaf) {
          std::cerr << "Removing: " << gene_ << " " << ts << " " << positions_[ts][i] << " | frequency: "
                    << maf[i] << " due to MAF filter." << std::endl;
        }
        sums.shed_col(i);
        genotypes_[ts].shed_col(i);
        missing_variant_carriers_[ts].shed_col(i);
        positions_[ts].erase(positions_[ts].begin() + i);
        nvariants_[ts]--;
      }
    }
    weights_[ts] = arma::vec(nvariants_[ts], arma::fill::zeros);
    // Check if any polymorphic. Mark transcripts skippable if all fixed.
    if (arma::accu(sums) == 0 || nvariants_[ts] == 0) {
      std::cerr << "gene: " << gene_ << " marked skippable." << std::endl;
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
      GLM<Gaussian> fit(D, Y, link);

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
      detail << gene_ << "\t";
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
      detail << gene_ << "\t";
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

auto Gene::get_detail() -> std::string {
  return detail_;
}

auto Gene::get_vaast() -> std::map<std::string, std::string> {
  return vaast_;
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
  if (nmaj > 0) {
    // There are more cases to distribute than minor allele carriers
    extreme_phen(maj_carriers(arma::span(0, nmaj - 1))).ones();
  }
  assert(arma::accu(extreme_phen) == ncase);

  VAAST vaast(genotypes_[k], extreme_phen, weights_[k], positions_[k], k, tp.score_only_minor, tp
      .score_only_alternative, false, tp.group_size, 2.);

  return arma::accu(vaast.expanded_scores > 0) >= 4;
}

auto Gene::is_testable() -> bool {
  return testable_;
}

auto Gene::is_skippable() -> bool {
  return skippable_;
}

auto Gene::is_polymorphic(const std::string &k) -> bool {
  return polymorphic_[k];
}

auto Gene::generate_vaast(Covariates &cov) -> void {
  // One entry per transcript
  for (const auto &ts : transcripts_) {
    std::stringstream vaast_ss;
    vaast_ss << ">" << ts << "\t" << gene_ << std::endl;

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
      vaast_ss << "N|N\t";

      if (case_alt > 0 && cont_alt == 0) {
        vaast_ss << "N|";
      } else if (case_alt > 0 && cont_alt > 0) {
        vaast_ss << "B|";
      }
      arma::uvec carriers = arma::find(X.col(i) % Y > 0);

      for (arma::uword j = 0; j < carriers.n_elem; j++) {
        if (j < carriers.n_elem - 1) {
          vaast_ss << carriers(j) << ",";
        } else {
          vaast_ss << carriers(j) << "|";
        }
      }
      // Placeholder
      vaast_ss << "N:N|N:N" << std::endl;
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
      vaast_ss << "N|N\t";

      arma::uvec carriers = arma::find(X.col(i) % (1. - Y) > 0);

      for (arma::uword j = 0; j < carriers.n_elem; j++) {
        if (j < carriers.n_elem - 1) {
          vaast_ss << carriers(j) << ",";
        } else {
          vaast_ss << carriers(j) << "|";
        }
      }
      // Placeholder
      vaast_ss << "N:N|N:N" << std::endl;
    }
    vaast_[ts] = vaast_ss.str();
  }
}

arma::sp_mat &Gene::get_missing(const std::string &k) {
  return missing_variant_carriers_[k];
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

