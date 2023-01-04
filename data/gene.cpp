//
// Created by Bohlender,Ryan James on 8/21/18.
//

#include <string>
#include <unordered_map>

#include <boost/format.hpp>
#include <stack>

#include "../link/binomial.hpp"
#include "../statistics/glm.hpp"
#include "../statistics/methods.hpp"
#include "../utility/math.hpp"
#include "filter.hpp"
#include "gene.hpp"

Gene::Gene(std::stringstream &ss, const std::shared_ptr<Covariates> &cov,
           unsigned long nsamples,
           std::unordered_map<std::string, arma::uword> nvariants,
           const Weight &weight, const TaskParams &tp, Filter &filter)
    : nsamples(nsamples), nvariants(std::move(nvariants)), testable(false),
      skippable(false), tp(tp) {
  parse(ss, cov, filter);

  if (tp.impute_to_mean) {
    impute_to_mean(cov);
  }

  if (tp.aaf_filter) {
    aaf_filter();
  } else {
    maf_filter();
  }

  update_weights(weight);
}

void Gene::aaf_filter() { // Mark gene skippable to start
  skippable = true;
  // Drop variants with minor allele count > tp_.mac.
  for (const auto &ts : transcripts) {
    arma::rowvec sums = arma::rowvec(arma::sum(genotypes[ts], 0));

    arma::rowvec aaf = arma::rowvec(arma::mean(genotypes[ts]) / 2.);
    for (arma::sword k = sums.n_elem - 1; k >= 0; --k) {
      bool bmac = sums[k] > tp.mac;
      bool baaf = aaf[k] >= .5;
      bool bwht = false; // whitelist bool
      if (!to_remove[ts].empty()) {
        bwht = k == to_remove[ts].top();
      }
      if (bmac || baaf || bwht) {
        if (tp.verbose && bmac) {
          std::cerr << "Removing: " << gene_name << " " << ts << " "
                    << positions[ts][k] << " | count: " << sums[k]
                    << " due to MAC filter." << std::endl;
        } else if (tp.verbose && baaf) {
          std::cerr << "Removing: " << gene_name << " " << ts << " "
                    << positions[ts][k] << " | frequency: " << aaf[k];
          if (tp.aaf_filter) {
            std::cerr << " due to AAF filter." << std::endl;
          } else {
            std::cerr << " due to MAF filter." << std::endl;
          }
        } else if (tp.verbose && bwht) {
          std::cerr << "Removing: " << gene_name << " " << ts << " "
                    << positions[ts][k] << " | type: " << type[ts][k]
                    << " | function: " << function[ts][k]
                    << " due to variant whitelist." << std::endl;
        }
        if (bwht) {
          to_remove[ts].pop();
        }
        sums.shed_col(k);
        genotypes[ts].shed_col(k);
        weights[ts].shed_row(k);
        missing_variant_carriers_[ts].shed_col(k);
        positions[ts].erase(positions[ts].begin() + k);
        function[ts].erase(function[ts].begin() + k);
        reference[ts].erase(reference[ts].begin() + k);
        alternate[ts].erase(alternate[ts].begin() + k);
        annotation[ts].erase(annotation[ts].begin() + k);
        type[ts].erase(type[ts].begin() + k);
        nvariants[ts]--;
      }
    }
    // Check if any polymorphic. Mark transcripts skippable if all fixed.
    if (arma::accu(sums) < tp.min_minor_allele_count ||
        nvariants[ts] < tp.min_variant_count) {
      std::cerr << "transcript: " << ts << " marked skippable." << std::endl;
      polymorphic[ts] = false;
    } else {
      polymorphic[ts] = true;
      skippable = false;
    }
  }
}

void Gene::maf_filter() {
  // Switch to counting minor allele
  for (auto &v : genotypes) {
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
  skippable = true;
  // Drop variants with minor allele count > tp_.mac.
  for (const auto &ts : transcripts) {
    arma::rowvec sums = arma::rowvec(arma::sum(genotypes[ts], 0));

    arma::rowvec maf = arma::rowvec(arma::mean(genotypes[ts]) / 2.);
    for (arma::sword k = sums.n_elem - 1; k >= 0; --k) {
      bool bmac = sums[k] > tp.mac;
      bool bmaf = maf[k] > tp.maf;
      bool bwht = false; // whitelist bool
      if (!to_remove[ts].empty()) {
        bwht = k == to_remove[ts].top();
      }
      if (bmac || bmaf || bwht) {
        if (tp.verbose && bmac) {
          std::cerr << "Removing: " << gene_name << " " << ts << " "
                    << positions[ts][k] << " | count: " << sums[k]
                    << " due to MAC filter." << std::endl;
        } else if (tp.verbose && bmaf) {
          std::cerr << "Removing: " << gene_name << " " << ts << " "
                    << positions[ts][k] << " | frequency: " << maf[k]
                    << " due to MAF filter." << std::endl;
        } else if (tp.verbose && bwht) {
          std::cerr << "Removing: " << gene_name << " " << ts << " "
                    << positions[ts][k] << " | type: " << type[ts][k]
                    << " | function: " << function[ts][k]
                    << " due to variant whitelist." << std::endl;
        }
        if (bwht) {
          to_remove[ts].pop();
        }
        sums.shed_col(k);
        genotypes[ts].shed_col(k);
        weights[ts].shed_row(k);
        missing_variant_carriers_[ts].shed_col(k);
        positions[ts].erase(positions[ts].begin() + k);
        function[ts].erase(function[ts].begin() + k);
        reference[ts].erase(reference[ts].begin() + k);
        alternate[ts].erase(alternate[ts].begin() + k);
        annotation[ts].erase(annotation[ts].begin() + k);
        type[ts].erase(type[ts].begin() + k);
        nvariants[ts]--;
      }
    }
    // Check if any polymorphic. Mark transcripts skippable if all fixed.
    if (arma::accu(sums) < tp.min_minor_allele_count ||
        nvariants[ts] < tp.min_variant_count) {
      std::cerr << "gene: " << gene_name << " marked skippable." << std::endl;
      polymorphic[ts] = false;
    } else {
      polymorphic[ts] = true;
      skippable = false;
    }
  }
}

void Gene::impute_to_mean(const std::shared_ptr<Covariates> &cov) {
  // Impute missing values to the mean of their category
  for (auto &[_, g] : genotypes) {
    arma::mat X(g);
    arma::vec Y(cov->get_original_phenotypes());
    arma::uvec cases = arma::find(Y == 1);
    arma::uvec controls = arma::find(Y == 0);
    for (arma::uword i = 0; i < X.n_cols; i++) {
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
    g = arma::sp_mat(X);
  }
}

void Gene::update_weights(const Weight &weight) {
  if (!weight.empty()) {
    for (const auto &k : transcripts) {
      weights[k].reshape(nvariants[k], 1);
      arma::uword i = 0;
      for (const auto &v : positions[k]) {
        weights[k](i) = weight.get(v);
        i++;
      }
    }
  } else {
    for (const auto &k : transcripts) {
      weights[k].reshape(nvariants[k], 1);
      arma::uword i = 0;
      for (const auto &v : positions[k]) {
        weights[k](i) = 1;
        i++;
      }
    }
  }
}

void Gene::print() {
  for (const auto &v : genotypes)
    std::cout << v.second;
}

void Gene::set_matrix(const std::string &k, arma::sp_mat &data) {
  genotypes[k] = std::move(data);
}

void Gene::set_matrix(const std::string &k, arma::sp_mat &&data) {
  genotypes[k] = std::move(data);
}

std::vector<std::string> &Gene::get_transcripts() { return transcripts; }

arma::vec &Gene::get_weights(const std::string &k) { return weights[k]; }

std::vector<std::string> &Gene::get_positions(const std::string &k) {
  return positions[k];
}

void Gene::clear(Covariates &cov,
                 std::unordered_map<std::string, Result> &results,
                 const TaskParams &tp) {
  if (!tp.no_detail)
    generate_detail(cov, results);
  if (tp.method == "VAAST") {
    generate_vaast(cov);
  }

  // Set matrix size to 0x0 to free space.
  for (auto &v : genotypes) {
    v.second.reset();
  }
}

void Gene::parse(std::stringstream &ss, const std::shared_ptr<Covariates> &cov,
                 Filter &filter) {
  using namespace RJBUtil;
  std::string line;
  arma::uword i = 0;
  std::vector<bool> rarest;

  if (tp.whole_gene) {
    ss = transcript_union(ss, cov, filter);
  }

  while (std::getline(ss, line, '\n')) {
    if (i == 0) {
      Splitter<std::string> splitter(line, "\t");
      for (arma::uword j = static_cast<int>(Indices::first);
           j < splitter.size(); j++) {
        if (cov->contains(splitter[j])) {
          samples.emplace_back(splitter[j]);
          columns.push_back(j);
        }
      }
      i++;
      continue;
    }
    if (columns.size() != cov->get_nsamples()) {
      std::cerr
          << "ERROR: Not all samples in .ped file are represented in the "
             "matrix file. There may be an ID mismatch or missing samples.";
      std::exit(-1);
    }
    Splitter<std::string> splitter(line, "\t", 11);

    if (gene_name.empty()) {
      gene_name = splitter[static_cast<int>(Indices::gene)];
    }
    std::string transcript =
        tp.whole_gene
            ? gene_name
            : std::string(splitter[static_cast<int>(Indices::transcript)]);
    auto iter =
        std::find(std::begin(transcripts), std::end(transcripts), transcript);
    if (iter == std::end(transcripts)) {
      // Transcript not found -- add
      transcripts.push_back(transcript);
      // Start with matrix transposed
      genotypes.emplace(std::make_pair(
          transcript, arma::sp_mat(nsamples, nvariants[transcript])));
      // Separately record the controls
      missing_variant_carriers_.emplace(std::make_pair(
          transcript, arma::sp_mat(nsamples, nvariants[transcript])));
      // Ensure positions are allocated for current transcript
      positions[transcript] = std::vector<std::string>();
      // Ensure weights are the right length
      weights[transcript].reshape(nvariants[transcript], 1);
      // Add transcript to the type map
      type.emplace(std::make_pair(transcript, std::vector<std::string>()));
      // Add transcript to the reference allele map
      reference.emplace(std::make_pair(transcript, std::vector<std::string>()));
      // Add transcript to the alternate allele map
      alternate.emplace(std::make_pair(transcript, std::vector<std::string>()));
      // Add transcript to the variant function map
      function.emplace(std::make_pair(transcript, std::vector<std::string>()));
      // Add transcript to the variant annotation map
      annotation.emplace(
          std::make_pair(transcript, std::vector<std::string>()));
      // Add filtering stack for transcript
      to_remove.emplace(std::make_pair(transcript, std::stack<int>()));
      // Reset counter on new transcript
      i = 1;
    }
    if (!filter.allow_variant(tp.method, splitter)) {
      to_remove[transcript].push(i - 1);
    }

    std::string var_id = form_variant_id(splitter);
    positions[transcript].push_back(var_id);

    // Record other info in line
    function[transcript].emplace_back(
        splitter[static_cast<int>(Indices::function)]);
    annotation[transcript].emplace_back(
        splitter[static_cast<int>(Indices::annotation)]);
    reference[transcript].emplace_back(
        splitter[static_cast<int>(Indices::ref)]);
    alternate[transcript].emplace_back(
        splitter[static_cast<int>(Indices::alt)]);
    type[transcript].emplace_back(splitter[static_cast<int>(Indices::type)]);

    auto first_idx = static_cast<int>(Indices::first);
    arma::uword k = 0;
    for (const auto &j : columns) {
      double val;
      try {
        val = splitter.back()[j - static_cast<int>(Indices::first)] - '0';
      } catch (std::exception &e) {
        std::cerr << "Full buffer: " << ss.str() << std::endl;
        std::cerr << "Failed to convert data to double: " << splitter[j]
                  << std::endl;
        std::cerr << "Line: " << line << std::endl;
        std::cerr << "j: " << j << std::endl;
        std::exit(-1);
      }
      if (val == -48) {
        std::cerr << splitter.back()[j - static_cast<int>(Indices::first)] << std::endl;
      }
      // Handle missing data
      if (val > 2) {
        if (tp.impute_to_mean) {
          val = -9;
        } else {
          val = 0;
          missing_variant_carriers_[transcript](k, i - 1) = 1;
        }
      }
      if (val != 0) {
        try {
          genotypes[transcript](k, i - 1) = val;
        } catch (std::exception &e) {
          std::cerr << "Failed to set value for: row = " << j - first_idx
                    << " col = " << i - 1 << std::endl;
          std::cerr << "Gene: " << splitter[static_cast<int>(Indices::gene)]
                    << " Transcript: "
                    << splitter[static_cast<int>(Indices::transcript)]
                    << " Location: " << form_variant_id(splitter) << std::endl;
          std::cerr << e.what();
          throw(e);
        }
      }
      k++;
    }
    i++;
  }
}

void Gene::set_weights(const std::string &k, arma::vec &new_weights) {
  if (new_weights.n_rows != genotypes[k].n_cols) {
    throw(std::logic_error("Weights do not match number of variants."));
  }

  weights[k] = new_weights;
}

void Gene::set_scores(const std::string &k, arma::vec &scores) {
  variant_scores[k] = scores;
}

void Gene::generate_detail(Covariates &cov,
                           std::unordered_map<std::string, Result> &results) {
  std::stringstream detail;

  arma::vec Y = cov.get_original_phenotypes();

  std::unordered_map<std::string, std::vector<std::string>> pos_ts_map;
  std::unordered_map<std::string, double> pos_score_map;
  std::unordered_map<std::string, double> pos_weight_map;
  std::unordered_map<std::string, double> pos_freq_map;

  // Case/Control Ref and Alt counts
  std::unordered_map<std::string, double> pos_caseref_map;
  std::unordered_map<std::string, double> pos_casealt_map;
  std::unordered_map<std::string, double> pos_contref_map;
  std::unordered_map<std::string, double> pos_contalt_map;

  // Alt carriers
  std::unordered_map<std::string, arma::uvec> pos_caseidx_map;
  std::unordered_map<std::string, arma::uvec> pos_contidx_map;

  // Collect all positions across transcripts and associated scores
  for (const auto &ts : transcripts) {
    arma::uword i = 0;

    arma::uvec cases = arma::find(Y == 1);
    arma::uvec controls = arma::find(Y == 0);

    arma::sp_mat X(genotypes[ts]);

    arma::rowvec maf = arma::rowvec(arma::mean(X) / 2.);
    // Ref/Alt Counts
    arma::vec case_alt = X.t() * Y;
    arma::vec case_ref = 2 * cases.n_elem - case_alt;
    arma::vec cont_alt = X.t() * (1. - Y);
    arma::vec cont_ref = 2 * controls.n_elem - cont_alt;
    if (!tp.qtl) {
      for (const auto &pos : positions[ts]) {
        // Get transcripts
        if (pos_ts_map.find(pos) == pos_ts_map.end()) {
          pos_ts_map[pos] = {ts};
        } else {
          pos_ts_map[pos].push_back(ts);
        }
        // Get scores
        if (pos_score_map.find(pos) == pos_score_map.end()) {
          if (variant_scores[ts].empty())
            variant_scores[ts].zeros(positions[ts].size());
          pos_score_map[pos] = variant_scores[ts](i);
          pos_weight_map[pos] = weights[ts](i);
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
        results[ts].testable = check_testability(ts, cov, results[ts].permuted);
      if (!testable)
        testable = results[ts].testable;
    } else {
      for (const auto &pos : positions[ts]) {
        // Get transcripts
        if (pos_ts_map.find(pos) == pos_ts_map.end()) {
          pos_ts_map[pos] = {ts};
        } else {
          pos_ts_map[pos].push_back(ts);
        }
        // Get scores
        if (pos_score_map.find(pos) == pos_score_map.end()) {
          if (variant_scores[ts].empty())
            variant_scores[ts].zeros(positions[ts].size());
          pos_score_map[pos] = variant_scores[ts](i);
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
  if (tp.qtl) {
    // Output for quantitative traits
    for (const auto &v : pos_ts_map) {
      detail << gene_name << "\t";
      print_comma_sep(v.second, detail);
      detail << "\t";
      detail << boost::format("%1$s\t%2$.2f\t%3$.2f\t%4$d") % v.first %
                    pos_score_map[v.first] % pos_weight_map[v.first] %
                    pos_freq_map[v.first];
      detail << std::endl;
    }
  } else {
    // Output for binary traits
    for (const auto &v : pos_ts_map) {
      detail << gene_name << "\t";
      print_comma_sep(v.second, detail);
      detail << "\t";
      detail << boost::format(
                    "%1$s\t%2$.2f\t%3$.2f\t%4$d\t%5$d\t%6$d\t%7$d\t%8$d") %
                    v.first % pos_score_map[v.first] % pos_weight_map[v.first] %
                    pos_freq_map[v.first] % pos_caseref_map[v.first] %
                    pos_casealt_map[v.first] % pos_contref_map[v.first] %
                    pos_contalt_map[v.first];
      detail << "\t";
      print_semicolon_sep(pos_caseidx_map[v.first], detail);
      detail << "\t";
      print_semicolon_sep(pos_contidx_map[v.first], detail);
      detail << std::endl;
    }
  }
  detail_ = detail.str();
}

std::string Gene::get_detail() const { return detail_; }

std::unordered_map<std::string, std::string> Gene::get_vaast() const {
  return vaast_;
}

const std::vector<std::string> &Gene::get_samples() const { return samples; }

bool Gene::check_testability(const std::string &transcript, Covariates &cov,
                             const std::vector<double> &permuted) {
  if (!is_polymorphic(transcript)) {
    return false;
  }
  double alpha = *tp.testable;
  arma::vec Yvec = cov.get_original_phenotypes();
  arma::sp_mat Xmat(genotypes[transcript]);

  auto ncase = static_cast<arma::sword>(arma::accu(Yvec));

  Methods methods(tp, cov);

  // Most Extreme Phenotype Distribution
  arma::uvec mac_carriers = arma::vec(arma::sum(Xmat, 1)) > 0;
  arma::uvec maj_carriers = 1 - mac_carriers;

  arma::sword nmac =
      std::min(ncase, static_cast<arma::sword>(arma::accu(mac_carriers)));
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

  Gene tmp = *this;
  double score = methods.call(tmp, cov, extreme_phen, transcript, true);

  if (tp.method == "VAAST" || tp.method == "SKAT") {
    bool testable_variants =
        arma::accu(tmp.variant_scores[transcript] > 0) >= 8;
    if (tp.analytic) {
      return score < alpha && testable_variants;
    } else {
      return (percentile_of_score<double>(score, permuted, false) < alpha) &&
             testable_variants;
    }
  } else {
    if (tp.analytic) {
      return score < alpha;
    } else {
      return percentile_of_score(score, permuted, false) < alpha;
    }
  }
}

bool Gene::is_skippable() const { return skippable; }

bool Gene::is_polymorphic(const std::string &k) { return polymorphic[k]; }

void Gene::generate_vaast(Covariates &cov) {
  // One entry per transcript
  for (const auto &ts : transcripts) {
    std::stringstream vaast_ss;
    vaast_ss << ">" << ts << "\t" << gene_name << std::endl;

    arma::sp_mat X(genotypes[ts]);
    arma::vec Y(cov.get_original_phenotypes());

    arma::uvec cases = arma::find(Y == 1);
    arma::uvec controls = arma::find(Y == 0);
    // Variant Format:
    // Variants are nucleotide then amino acid
    // R/U: variant_score position@chr nt|AA B/N|id1,id2...|ref:alt|ref:alt

    // Skip exon positions for now
    // /Volumes/huff/yyu4/Analysis_case_control/caper/LUNG/lung.bing_0.0001
    // lung_cov_sex_age_pc
    // /rsrch3/home/epi/yyu4/workspace5/Analysis_case_control/caper/LUNG/lung_no_cov/c_18.b_1/VAAST.g2.vaast

    // Add target variants
    for (int i = 0; i < positions[ts].size(); i++) {
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
      vaast_ss << boost::format("%1$.2f") % variant_scores[ts][i] << "\t";

      RJBUtil::Splitter<std::string> pos_splitter(positions[ts][i], ",");
      vaast_ss << pos_splitter[1] << "@" << pos_splitter[0] << "\t";

      // Placeholder for Reference alleles
      vaast_ss << boost::format("%1$s\t") % reference[ts][i];
      vaast_ss << boost::format("%1$s\t") % annotation[ts][i];

      arma::vec Xcol(X.col(i));

      arma::uvec het_carriers = arma::find(Xcol % Y == 1);
      arma::uvec hom_carriers = arma::find(Xcol % Y == 2);

      if (het_carriers.n_elem > 0) {
        vaast_ss << compress_adjacent(het_carriers);
        // Variant
        vaast_ss << boost::format("%1$s:%2$s") % reference[ts][i] %
                        alternate[ts][i]
                 << " ";
      }
      if (hom_carriers.n_elem > 0) {
        vaast_ss << "\t";
        vaast_ss << compress_adjacent(hom_carriers);
        // Variant
        vaast_ss << boost::format("%1$s:%2$s") % alternate[ts][i] %
                        alternate[ts][i];
      }
      vaast_ss << std::endl;
    }
    // Add background variants
    for (int i = 0; i < positions[ts].size(); i++) {
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
      vaast_ss << boost::format("%1$.2f") % variant_scores[ts][i] << "\t";

      RJBUtil::Splitter<std::string> pos_splitter(positions[ts][i], ",");
      vaast_ss << pos_splitter[1] << "@" << pos_splitter[0] << "\t";

      // Placeholder for Reference alleles
      vaast_ss << boost::format("%1$s\t") % reference[ts][i];
      vaast_ss << boost::format("%1$s\t") % annotation[ts][i];

      arma::vec Xcol(X.col(i));

      arma::uvec het_carriers = arma::find(Xcol % (1. - Y) == 1);
      arma::uvec hom_carriers = arma::find(Xcol % (1. - Y) == 2);

      if (het_carriers.n_elem > 0) {
        for (arma::uword j = 0; j < het_carriers.n_elem; j++) {
          if (j < het_carriers.n_elem - 1) {
            vaast_ss << het_carriers(j) << ",";
          } else {
            vaast_ss << het_carriers(j) << "|";
          }
        }
        // Variant
        vaast_ss << boost::format("%1$s:%2$s") % reference[ts][i] %
                        alternate[ts][i]
                 << " ";
      }
      if (hom_carriers.n_elem > 0) {
        vaast_ss << "\t";
        for (arma::uword j = 0; j < hom_carriers.n_elem; j++) {
          if (j < hom_carriers.n_elem - 1) {
            vaast_ss << hom_carriers(j) << ",";
          } else {
            vaast_ss << hom_carriers(j) << "|";
          }
        }
        // Variant
        vaast_ss << boost::format("%1$s:%2$s") % alternate[ts][i] %
                        alternate[ts][i];
      }
      vaast_ss << std::endl;
    }
    vaast_[ts] = vaast_ss.str();
  }
}

std::string Gene::form_variant_id(RJBUtil::Splitter<std::string> &splitter) {
  std::stringstream ss;
  ss << splitter[static_cast<int>(Indices::chrom)] << ","
     << splitter[static_cast<int>(Indices::start)] << ","
     << splitter[static_cast<int>(Indices::end)] << ","
     << splitter[static_cast<int>(Indices::ref)] << ","
     << splitter[static_cast<int>(Indices::alt)] << ","
     << splitter[static_cast<int>(Indices::type)] << ","
     << splitter[static_cast<int>(Indices::gene)] << ","
     << splitter[static_cast<int>(Indices::transcript)];
  return ss.str();
}

std::stringstream Gene::transcript_union(std::stringstream &ss,
                                         const std::shared_ptr<Covariates> &cov,
                                         Filter &filter) {
  std::map<std::string, std::string> variants;
  std::stringstream res_ss;
  std::string line;
  arma::uword line_no = 0;

  std::map<std::string, double> casm_scores;

  while (std::getline(ss, line, '\n')) {
    if (line_no == 0) {
      line_no++;
      res_ss << line << "\n";
      continue;
    }
    RJBUtil::Splitter<std::string> splitter(line, "\t");
    std::string vid = form_variant_id(splitter);
    if (gene_name.empty()) {
      gene_name = splitter[static_cast<int>(Indices::gene)];
    }
    if (variants.find(vid) == variants.end()) {
      variants[vid] = line;
    }
    line_no++;
  }

  nvariants[gene_name] = variants.size();

  for (const auto &s : variants) {
    res_ss << s.second << "\n";
  }

  return res_ss;
}

std::string Gene::compress_adjacent(arma::uvec &samples) {
  std::stringstream ss;
  int i = 0;
  int j = 1;
  while (i < samples.n_elem) {
    int n = 1;
    j = i + 1;
    while (j < samples.n_elem && samples(i) + n == samples(j)) {
      j++;
      n++;
    }
    if (n > 1) {
      if (j < samples.n_elem - 1) {
        ss << samples(i) << "-" << samples(j - 1) << ",";
      } else {
        ss << samples(i) << "-" << samples(j - 1) << "|";
      }
    } else {
      if (i < samples.n_elem - 1) {
        ss << samples(i) << ",";
      } else {
        ss << samples(i) << "|";
      }
    }
    i = j;
  }
  return ss.str();
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
