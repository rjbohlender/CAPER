//
// Created by Bohlender,Ryan James on 7/31/18.
//

#define ARMA_DONT_PRINT_ERRORS

#include <cmath>
#include <iomanip>

// Boost Math
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/constants/constants.hpp>

#include <boost/format.hpp>

#include "methods.hpp"
#include "vaast.hpp"

#include "../link/binomial.hpp"
#include "../link/gaussian.hpp"
#include "fishertest.hpp"
#include "glm.hpp"

constexpr std::array<double, 8> Methods::rho_;

/**
 * @brief Compute the rank of each element in a vector
 * @param v The vector to rank
 * @param direction The direction to rank in. Either "ascend" or "descend"
 * @param method Tie breaking method, 0 = average, 1 = min, 2 = max
 * @return A vector of ranks
 */
arma::vec rank(arma::vec &v, const char *direction, int method) {
  if (strcmp(direction, "ascend") != 0 && strcmp(direction, "descend") != 0)
    throw(std::logic_error(
        "Order argument for rank() must be either 'ascend' or 'descend'"));

  arma::uvec sort_indices;
  try {
    sort_indices = arma::sort_index(v, direction);
  } catch (const std::logic_error &e) {
    std::cerr << "NANs among ranked values. Replacing with 0.\n";
    v.replace(arma::datum::nan, 0);
    sort_indices = arma::sort_index(v, direction);
  }
  arma::vec sorted = v(sort_indices);

  arma::vec ranks = arma::vec(v.n_rows, arma::fill::zeros);
  arma::sword i = 0, j = 0;

  while (i < v.n_rows) {
    j = i + 1;
    // Find the next different value
    while (j < v.n_rows) {
      if (sorted(i) != sorted(j))
        break;
      j++;
    }
    // Adjusted rank
    if (method == 0) { // average
      for (arma::uword k = i; k < j; k++) {
        ranks(sort_indices(k)) = 1. + (i + j - 1.) / 2.0f;
      }
    } else if(method == 1) { // min
      for (arma::uword k = i; k < j; k++) {
        ranks(sort_indices(k)) = i;
      }

    } else if(method == 2) { // max
      for (arma::uword k = i; k < j; k++) {
        ranks(sort_indices(k)) = j;
      }
    }
    // Update i
    i = j;
  }

  return ranks;
}

Methods::Methods(const TaskParams &tp, const std::shared_ptr<Covariates> &cov)
    : method_(tp.method), tp(tp) {
  if (tp.kernel == "Linear") {
    kernel_ = Kernel::Linear;
  } else if (tp.kernel == "wLinear") {
    kernel_ = Kernel::wLinear;
  }

  if ((tp.method == "SKAT" || tp.method == "SKATO" || tp.method == "BURDEN" || tp.method == "SKATC") &&
      !tp.qtl) {
    obj_ = std::make_shared<SKATR_Null>(*cov);
    lin_obj_ = nullptr;
  } else if ((tp.method == "SKAT" || tp.method == "SKATO" ||
              tp.method == "BURDEN" || tp.method == "SKATC") &&
             tp.qtl) {
    obj_ = nullptr;
    lin_obj_ = std::make_shared<SKATR_Linear_Null>(*cov);
  } else if (tp.method == "VT") {
    vt_obj_ = std::make_shared<VT_Res>();
  }
}

Methods::Methods(const TaskParams &tp, const Covariates &cov)
    : method_(tp.method), tp(tp) {
  if (tp.kernel == "Linear") {
    kernel_ = Kernel::Linear;
  } else if (tp.kernel == "wLinear") {
    kernel_ = Kernel::wLinear;
  }

  if ((tp.method == "SKAT" || tp.method == "SKATO" || tp.method == "BURDEN" || tp.method == "SKATC") &&
      !tp.qtl) {
    obj_ = std::make_shared<SKATR_Null>(cov);
    lin_obj_ = nullptr;
  } else if ((tp.method == "SKAT" || tp.method == "SKATO" ||
              tp.method == "BURDEN" || tp.method == "SKATC") &&
             tp.qtl) {
    obj_ = nullptr;
    lin_obj_ = std::make_shared<SKATR_Linear_Null>(cov);
  } else if (tp.method == "VT") {
    vt_obj_ = std::make_shared<VT_Res>();
  }
}

/**
 * @brief Reset kernel to free memory.
 */
void Methods::clear(std::vector<std::string> &v) {
  obj_.reset();
  lin_obj_.reset();
}

/**
 * @brief Compute the p-value using the BURDEN test statistic.
 * @param gene The gene to test.
 * @param ts Transcript.
 * @param phenotypes Phenoypes.
 * @param permute Whether we are permuting or returning analytic p-value.
 * @return The p-value or test statistic.
 */
double Methods::BURDEN(Gene &gene, arma::vec &phenotypes, const std::string &ts,
                       bool permute) {
  obj_->shuffle(phenotypes);
  arma::sp_mat G(gene.genotypes[ts]);
  arma::uword N = G.n_cols; // Variant count

  check_weights(gene, ts, tp.a, tp.b, tp.no_weights);
  arma::vec weights = gene.get_weights(ts);

  arma::mat W = arma::diagmat(weights);

  arma::mat tmp;
  if (tp.qtl) {
    tmp = lin_obj_->get_Ux().t() * G;
  } else {
    tmp = obj_->get_Ux().t() * G;
  }

  arma::mat Gs;
  arma::rowvec Zs;
  if (tp.qtl) {
    Gs = G.t() * G - tmp.t() * tmp;
    Zs = arma::sum(arma::diagmat(lin_obj_->get_U0()) * G) /
         std::sqrt(lin_obj_->get_s2());
  } else {
    Gs = (arma::diagmat(obj_->get_Yv()) * G).t() * G - tmp.t() * tmp;
    Zs = arma::sum(arma::diagmat(obj_->get_U0()) * G);
  }

  arma::mat R = (Gs * W).t() * W;
  arma::mat Z = Zs * W;

  double Qb = std::pow(arma::accu(Z), 2);

  arma::vec Rs = arma::sum(R, 1);
  double R1 = arma::accu(Rs);
  double stat = Qb / R1;

  if (!permute) {
    boost::math::chi_squared chisq(1); // 1-df chisq
    if (!std::isfinite(stat)) {
      return 1.;
    } else {
      return boost::math::cdf(boost::math::complement(chisq, stat));
    }
  } else {
    return stat;
  }
}

/**
 * @brief Calculate the p-value for a gene using the SKAT method.
 *
 * @param gene Gene object.
 * @param Y Phenotype vector.
 * @param ts Transcript.
 * @return double P-value.
 */
double Methods::CALPHA(Gene &gene, arma::vec &Y, const std::string &ts) {
  arma::sp_mat X(gene.genotypes[ts]);

  double nA = arma::sum(Y);  // Case count
  double nU = Y.n_rows - nA; // Control count

  double p0 = nA / (nA + nU);

  arma::vec n(X.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i < X.n_cols; i++) {
    for (const auto &v : X.col(i)) {
      if (v > 0)
        n(i)++;
    }
  }

  arma::vec g(X.n_cols, arma::fill::zeros);
  arma::uvec case_idx = arma::find(Y == 1);
  for (auto it = X.begin(); it != X.end(); ++it) {
    if (arma::find(case_idx == it.row()).eval().n_elem > 0) {
      g(it.col())++;
    }
  }
  // arma::vec g = arma::sum(X.rows(arma::find(Y == 1)) > 0, 0).t();

  // Test statistic
  return arma::sum(arma::pow(g - (n * p0), 2) - (n * p0 * (1 - p0)));
}

/**
 * @brief Compute the CMC test statistic.
 *
 * @param gene Gene object.
 * @param Y Phenotype vector.
 * @param ts Transcript
 * @param rare_freq Rare variant frequency threshold.
 * @return double Test statistic.
 */
double Methods::CMC(Gene &gene, arma::vec &Y, const std::string &ts,
                    double rare_freq) const {
  if (tp.hotellings) {
    arma::mat X(gene.genotypes[ts]);

    double N = Y.n_rows;
    double nA = arma::sum(Y); // Case count
    double nU = N - nA;       // Control count

    arma::rowvec MAF = arma::mean(X, 0) / 2;

    // Collapse rare variants
    arma::uvec rare = arma::find(MAF < rare_freq);
    arma::uvec common = arma::find(MAF >= rare_freq);

    arma::mat Xnew;
    if (rare.size() <= 1) {
      Xnew = X;
    } else {
      arma::mat Xcollapse = arma::sum(X.cols(rare), 1);
      Xcollapse(arma::find(Xcollapse > 1)).ones();
      Xnew = X.cols(common);
      Xnew.insert_cols(Xnew.n_cols, Xcollapse);
    }

    // Rescale to -1, 0, 1
    Xnew -= 1;

    // Calculate two-sample Hotelling's T2 statistic
    arma::mat Xx = Xnew.rows(arma::find(Y == 1));
    arma::mat Yy = Xnew.rows(arma::find(Y == 0));

    arma::rowvec Xxmean = arma::mean(Xx);
    arma::rowvec Yymean = arma::mean(Yy);

    arma::mat COV =
        ((nA - 1.) * arma::cov(Xx) + (nU - 1.) * arma::cov(Yy)) / (N - 2.);
    arma::mat INV;
    if (!arma::inv_sympd(INV, COV)) {
      arma::pinv(INV, COV);
    }
    double ret = arma::as_scalar((Xxmean - Yymean) * INV *
                                 (Xxmean - Yymean).t() * nA * nU / N);
    auto p = static_cast<double>(Xxmean.n_elem);
    double stat = ret * (nA + nU - 1 - p) / (p * (nA + nU - 2));
    // stat ~ F(N, nA + nU - 1 - N)
    if (stat < 0)
      stat = 0;
    if (tp.nperm > 0) {
      return stat;
    }
    boost::math::fisher_f fisher_f(p, nA + nU - 1 - p);
    double pval;
    try {
      if (isnan(stat)) {
        return 1.;
      }
      pval = boost::math::cdf(boost::math::complement(fisher_f, stat));
    } catch (boost::exception &e) {
      std::cerr << "COV: " << COV;
      std::cerr << "INV: " << INV;
      std::cerr << "Xxmean: " << Xxmean;
      std::cerr << "Yymean: " << Yymean;
      std::cerr << "stat: " << stat << std::endl;
      std::cerr << "ret: " << arma::as_scalar(ret);
      throw;
    }
    return pval;
  } else {
    arma::mat X(gene.genotypes[ts]);

    const arma::rowvec freq = arma::mean(X, 0) / 2.;

    // Collapse rare variants
    const arma::uvec rare = arma::find(freq < rare_freq);
    const arma::uvec common = arma::find(freq >= rare_freq);

    arma::mat Xnew;
    double ncollapse = rare.n_elem;
    if (ncollapse > 0) {
      const double combined_freq = arma::accu(freq(rare));
      arma::mat Xcollapse = arma::sum(X.cols(rare), 1);
      Xnew = X.cols(common);
      Xnew.insert_cols(Xnew.n_cols, Xcollapse);
    } else {
      Xnew = X;
    }

    if (Xnew.n_cols == 1) { // Switch test
      Xnew(arma::find(Xnew > 0)).ones();

      double freq_mutated = arma::accu(arma::mean(Xnew, 0));

      arma::uword ncase = arma::accu(Y);
      arma::uword ncont = arma::accu(1 - Y);

      // Observed
      double case_alt = arma::accu(Xnew % Y);
      double cont_alt = arma::accu(Xnew % (1 - Y));
      double case_ref = ncase - case_alt;
      double cont_ref = ncont - cont_alt;

      // Expected
      double case_alt_exp = ncase * freq_mutated;
      double case_ref_exp = ncase * (1. - freq_mutated);
      double cont_alt_exp = ncont * freq_mutated;
      double cont_ref_exp = ncont * (1. - freq_mutated);

      double stat = std::pow(case_alt - case_alt_exp, 2) / case_alt_exp +
                    std::pow(case_ref - case_ref_exp, 2) / case_ref_exp +
                    std::pow(cont_alt - cont_alt_exp, 2) / cont_alt_exp +
                    std::pow(cont_ref - cont_ref_exp, 2) / cont_ref_exp;
      double df = 1;

      if (tp.nperm == 0) {
        // Formula is df = (r -1)(c - 1) and r == 2 always, so it's just c - 1.
        boost::math::chi_squared chisq(df);
        if (stat < 0) {
          stat = std::numeric_limits<double>::epsilon();
        }
        return boost::math::cdf(boost::math::complement(chisq, stat));
      } else {
        return stat;
      }
    } else {

      double n = arma::accu(Xnew); // Total of table observations

      arma::rowvec case_counts = Y.t() * Xnew;
      arma::rowvec control_counts = (1. - Y).t() * Xnew;

      double case_total = arma::accu(case_counts);
      double control_total = arma::accu(control_counts);

      arma::rowvec variant_totals = case_counts + control_counts;

      arma::rowvec case_expected = case_total * variant_totals / n;
      arma::rowvec control_expected = control_total * variant_totals / n;

      arma::rowvec case_chi =
          arma::pow(case_counts - case_expected, 2) / case_expected;
      arma::rowvec control_chi =
          arma::pow(control_counts - control_expected, 2) / control_expected;

      if (tp.nperm == 0) {
        // Formula is df = (r -1)(c - 1) and r == 2 always, so it's just c - 1.
        const int df = case_counts.n_elem - 1;
        boost::math::chi_squared chisq(df);
        double stat = arma::accu(case_chi + control_chi);
        if (stat < 0) {
          stat = std::numeric_limits<double>::epsilon();
        }
        return boost::math::cdf(boost::math::complement(chisq, stat));
      } else {
        return arma::accu(case_chi + control_chi);
      }
    }
  }
}

/**
 * @brief Calculate the CMC1df statistic for a gene.
 *
 * @param gene Gene to calculate the statistic for.
 * @param Y Phenotype vector.
 * @param ts Transcript to calculate the statistic for.
 * @return double Statistic.
 */
double Methods::CMC1df(Gene &gene, arma::vec &Y, const std::string &ts,
                       bool permute) const {
  // Runtime for MDA OV with just fisher test and 10000 perms = 6544.95
  // Runtime for fast path with 10000 perms = 267.874
  if (permute) {
    arma::vec X(arma::sum(gene.genotypes[ts], 1));
    X(arma::find(X > 0)).ones();

    arma::uword ncase = arma::accu(Y);
    arma::uword ncont = arma::accu(1 - Y);

    double case_alt = arma::accu(X % Y);
    double cont_alt = arma::accu(X % (1 - Y));
    double case_ref = ncase - case_alt;
    double cont_ref = ncont - cont_alt;

    if (case_alt == 0 || cont_alt == 0 || case_ref == 0 || cont_ref == 0) {
      case_alt += 0.5;
      cont_alt += 0.5;
      case_ref += 0.5;
      cont_ref += 0.5;
    }

    double OR = case_alt * cont_ref / (cont_alt * case_ref);

    return OR;
  } else {
    FisherTest fisherTest(gene, Y, ts);
    return fisherTest.get_pval();
  }
}

double Methods::RVT1(Gene &gene, arma::vec &Y, arma::mat design,
                     arma::vec &initial_beta, const std::string &ts,
                     bool linear) {
  // Runtime 100 perms naive initialization on macbook pro, ovarian data --
  // 2421.09 Runtime 100 perms prior initialization on macbook pro, ovarian data
  // -- 1863.08 -- Poor initialization in permutation Runtime 100 perms single
  // prior initialization on macbook pro, ovarian data -- 1757.13
  double stat;
  if (linear) {
    // Quantitative trait
    const arma::sp_mat X = arma::ceil(gene.genotypes[ts].t() / 2);
    Gaussian link("identity");
    GLM<Gaussian> fit1(design, Y, link, initial_beta, tp);
    arma::mat d2 =
        arma::join_horiz(design, arma::rowvec(arma::sum(X) / X.n_rows).t());
    GLM<Gaussian> fit2(d2, Y, link, fit1.beta_, tp);

    const double n = Y.n_elem;
    if (tp.wald) {
      // Calculate the MSE of the fit
      const double mse = arma::accu(arma::pow(Y - fit2.mu_, 2)) / (n - d2.n_rows);
      arma::mat var_beta_ = mse * arma::inv(d2.t() * d2);
      const double var_beta = var_beta_(var_beta_.n_rows - 1, var_beta_.n_cols - 1);

      stat = fit2.beta_(fit2.beta_.n_elem - 1) / std::sqrt(var_beta);
    } else {
      stat = (fit1.dev_ - fit2.dev_) / (fit2.dev_ / (n - arma::rank(d2)));
    }
  } else {
    // Binary trait
    // Convert to 0/1 to make summing the number of carriers easier.
    arma::sp_mat X = arma::ceil(gene.genotypes[ts].t() / 2);
    Binomial link("logit");
    GLM<Binomial> fit1(design, Y, link, initial_beta, tp);
    arma::mat d2 =
        arma::join_horiz(design, arma::rowvec(arma::sum(X) / X.n_rows).t());
    GLM<Binomial> fit2(d2, Y, link, fit1.beta_, tp);

    if (tp.wald) {
      // Calculate the MSE of the fit
      arma::mat var_beta_ = arma::inv(d2.t() * arma::diagmat(fit2.mu_ % (1. - fit2.mu_)) * d2);
      const double var_beta = var_beta_(var_beta_.n_rows - 1, var_beta_.n_cols - 1);

      stat = std::pow(fit2.beta_(fit2.beta_.n_elem - 1), 2) / var_beta;
    } else {
      stat = fit1.dev_ - fit2.dev_;
    }
  }
  boost::math::chi_squared chisq(1);
  if (stat < 0) {
    stat = std::numeric_limits<double>::epsilon();
  }
  return boost::math::cdf(boost::math::complement(chisq, stat));
}

double Methods::RVT2(Gene &gene, arma::vec &Y, arma::mat design,
                     arma::vec &initial_beta, const std::string &ts,
                     bool linear) {
  double stat;
  if (linear) {
    // Quantitative trait
    const arma::sp_mat X = gene.genotypes[ts].t();
    Gaussian link("identity");
    GLM<Gaussian> fit1(design, Y, link, initial_beta, tp);
    arma::rowvec r =
        arma::conv_to<arma::rowvec>::from(arma::rowvec(arma::sum(X)) > 0);
    arma::mat d2 = arma::join_horiz(design, r.t());
    GLM<Gaussian> fit2(d2, Y, link, fit1.beta_, tp);

    double n = Y.n_elem;

    if (tp.wald) {
      // Calculate the MSE of the fit
      const double mse = arma::accu(arma::pow(Y - fit2.mu_, 2)) / (n - d2.n_rows);
      arma::mat var_beta_ = mse * arma::inv(d2.t() * d2);
      const double var_beta = var_beta_(var_beta_.n_rows - 1, var_beta_.n_cols - 1);

      stat = fit2.beta_(fit2.beta_.n_elem - 1) / std::sqrt(var_beta);
    } else {
      stat = (fit1.dev_ - fit2.dev_) / (fit2.dev_ / (n - arma::rank(d2)));
    }
  } else {
    // Binary trait
    arma::sp_mat X = gene.genotypes[ts].t();
    Binomial link("logit");
    GLM<Binomial> fit1(design, Y, link, initial_beta, tp);
    arma::rowvec r =
        arma::conv_to<arma::rowvec>::from(arma::rowvec(arma::sum(X)) > 0);
    arma::mat d2 = arma::join_horiz(design, r.t());
    GLM<Binomial> fit2(d2, Y, link, fit1.beta_, tp);

    if (tp.wald) {
      // Calculate the variance of the betas
      arma::mat var_beta_ = arma::inv(d2.t() * arma::diagmat(fit2.mu_ % (1. - fit2.mu_)) * d2);
      const double var_beta = var_beta_(var_beta_.n_rows - 1, var_beta_.n_cols - 1);

      stat = std::pow(fit2.beta_(fit2.beta_.n_elem - 1), 2) / var_beta;
    } else {
      stat = fit1.dev_ - fit2.dev_;
    }
  }
  boost::math::chi_squared chisq(1);
  if (stat < 0) {
    stat = std::numeric_limits<double>::epsilon();
  }
  return boost::math::cdf(boost::math::complement(chisq, stat));
}

double Methods::VAAST(Gene &gene, arma::vec &Y, const std::string &ts,
                      double site_penalty, arma::uword group_threshold,
                      bool detail, bool biallelic, double control_freq_cutoff,
                      bool legacy) {
  check_weights(gene, ts, tp.a, tp.b, tp.no_weights);
  VAASTLogic vaast_logic(gene, Y, ts, site_penalty, group_threshold, detail,
                         biallelic, control_freq_cutoff, legacy);
  return vaast_logic.score;
}

double Methods::VT(Gene &gene, arma::vec &phenotypes, const std::string &ts) {
  // Convert data to match their format
  arma::vec pheno = arma::repmat(phenotypes, gene.genotypes[ts].n_cols, 1);
  if (!vt_obj_->is_initialized(ts)) {
    vt_obj_->initialize(gene, pheno, ts);
  }
  arma::vec phenoCount =
      pheno % vt_obj_->get_mCount(ts); // Changes under permutation
  arma::vec csPhenoCount = arma::cumsum(vt_obj_->sum_groups(
      phenoCount, vt_obj_->get_oneToLen(ts), ts)); // Changes under permutation

  return arma::max((csPhenoCount - vt_obj_->get_csCountMeanpheno(ts)) /
                   vt_obj_->get_sqrtCsCountSquare(ts));
}

double Methods::WSS(Gene &gene, arma::vec &Y, const std::string &k) {
  arma::mat X(gene.genotypes[k]);

  double n = Y.n_rows;

  arma::vec count = arma::sum(X, 0).t();         // Total count for each variant
  arma::vec freq = (1. + count) / (2. + 2. * n); // Frequency with a prior
  arma::vec weight =
      1. / arma::sqrt(freq % (1. - freq)); // Reciprocal of Binomial SD
  arma::vec count_weight = arma::vec(X * weight);

  return arma::accu(count_weight % Y);
}

std::string Methods::str() { return method_; }

/**
 * @brief Calculate SKAT with p-value following Wu, Guan, Pankow (2017)
 * @param gene Gene object to calculate SKAT on
 * @param transcript Transcript to calculate SKAT on
 * @param phenotypes Phenotypes to calculate SKAT on
 * @param a Parameter for Beta distribution
 * @param b Parameter for Beta distribution
 * @param detail Whether to generate the full output
 * @param linear Whether to use linear regression
 * @param permute Whether to permute the phenotypes
 * @return p-value or test statistic
 */
double Methods::SKAT(Gene &gene, arma::vec &phenotypes,
                     const std::string &transcript, int a, int b, bool detail,
                     bool linear, bool permute) {
  arma::sp_mat G(gene.genotypes[transcript]);

  if (linear) {
    lin_obj_->shuffle(phenotypes);
  } else {
    obj_->shuffle(phenotypes);
  }

  check_weights(gene, transcript, a, b, tp.no_weights);
  arma::vec weights = gene.get_weights(transcript);

  arma::mat W = arma::diagmat(weights);
  // We're permuting, only calculate the Q-value
  if (permute) {
    arma::rowvec Zs;
    if (linear) {
      Zs = arma::sum(arma::diagmat(lin_obj_->get_U0()) * G) /
           std::sqrt(lin_obj_->get_s2());
    } else {
      Zs = arma::sum(arma::diagmat(obj_->get_U0()) * G);
    }
    arma::mat Z = Zs * W;

    double Q = arma::accu(arma::pow(Z, 2));

    return Q;
  } else {
    // We're not permuting, return asymptotic p-values
    arma::mat tmp;
    if (linear) {
      tmp = lin_obj_->get_Ux().t() * G;
    } else {
      tmp = obj_->get_Ux().t() * G;
    }

    arma::mat Gs;
    arma::rowvec Zs;
    if (linear) {
      Gs = G.t() * G - tmp.t() * tmp;
      Zs = arma::sum(arma::diagmat(lin_obj_->get_U0()) * G) /
           std::sqrt(lin_obj_->get_s2());
    } else {
      Gs = (arma::diagmat(obj_->get_Yv()) * G).t() * G - tmp.t() * tmp;
      Zs = arma::sum(arma::diagmat(obj_->get_U0()) * G);
    }

    arma::mat R = (Gs * W).t() * W;
    arma::mat Z = Zs * W;

    arma::vec s;
    arma::svd(s, R);

    double Q = arma::accu(arma::pow(Z, 2));

    if (detail) {
      arma::vec variant_scores = arma::sum(arma::pow(Z, 2), 0).t();
      gene.set_scores(transcript, variant_scores);
    }

    return SKAT_pval(Q, s, tp.saddlepoint);
  }
}

double Methods::SKATO(Gene &gene, arma::vec &phenotypes,
                      const std::string &transcript, int a, int b, bool detail,
                      bool linear) {
  if (linear) {
    lin_obj_->shuffle(phenotypes);
  } else {
    obj_->shuffle(phenotypes);
  }

  arma::mat Gtmp(gene.genotypes[transcript]);
  arma::uvec Gidx = arma::find(arma::sum(Gtmp, 0) < 10);
  // Collapse rare variants into a single variant
  arma::sp_mat G;
#if 1
  if (Gidx.n_elem > 0) {
    // Insert column containing the collapsed variant
    Gtmp.insert_cols(Gtmp.n_cols, arma::sum(Gtmp.cols(Gidx), 1));
    // Remove the rare variants
    Gtmp.shed_cols(Gidx);
    // Update the genotype matrix
    G = Gtmp;
  } else {
    G = Gtmp;
  }
#else
  G = Gtmp;
#endif
  arma::uword N = G.n_cols; // Variant count

  check_weights(gene, transcript, a, b, tp.no_weights);
  arma::vec weights = gene.get_weights(transcript);
#if 1
  weights.insert_rows(weights.n_elem, arma::sum(weights.rows(Gidx)));
  weights.shed_rows(Gidx);
#endif

  arma::mat W = arma::diagmat(weights);

  arma::mat tmp;
  if (linear) {
    tmp = lin_obj_->get_Ux().t() * G;
  } else {
    tmp = obj_->get_Ux().t() * G;
  }

  arma::mat Gs;
  arma::rowvec Zs;
  if (linear) {
    Gs = G.t() * G - tmp.t() * tmp;
    Zs = arma::sum(arma::diagmat(lin_obj_->get_U0()) * G) /
         std::sqrt(lin_obj_->get_s2());
  } else {
    Gs = (arma::diagmat(obj_->get_Yv()) * G).t() * G - tmp.t() * tmp;
    Zs = arma::sum(arma::diagmat(obj_->get_U0()) * G);
  }

  arma::mat R = (Gs * W).t() * W;
  arma::mat Z = Zs * W;

  arma::vec s;
  arma::svd(s, R);

  arma::uword K = 8; // Length of rho_

  double Qs = arma::accu(arma::pow(Z, 2));
  double Qb = std::pow(arma::accu(Z), 2);
  arma::vec Qw{0, 0, 0, 0, 0, 0, 0, 0};

  for (arma::uword i = 0; i < K; i++) {
    Qw[i] = (1 - rho_[i]) * Qs + rho_[i] * Qb;
  }

  arma::vec pval = {0, 0, 0, 0, 0, 0, 0, 0};

  arma::vec Rs = arma::sum(R, 1);
  double R1 = arma::accu(Rs);
  double R2 = arma::accu(arma::pow(Rs, 2));
  // double R3 = arma::accu(Rs % arma::sum(R.each_col() % Rs).t());
  double R3 = arma::accu(Rs.t() * R * Rs);

  arma::mat RJ2(Rs.n_rows, Rs.n_rows, arma::fill::zeros);
  for (arma::uword i = 0; i < Rs.n_rows; i++) {
    RJ2.row(i) = (Rs(i) + Rs.t()) / N;
  } // Replacement for R's outer(Rs, Rs, '+')

  std::vector<arma::vec> lamk(K - 1);
  for (arma::uword i = 0; i < K; i++) {
    // Pure burden
    if (rho_[i] == 1) {
      boost::math::chi_squared chisq(1); // 1-df chisq
      double stat = Qb / R1;
      if (!std::isfinite(stat)) {
        pval[i] = 1.;
      } else {
        pval[i] = boost::math::cdf(boost::math::complement(chisq, stat));
      }
      continue;
    }

    // Setup Davies
    double c1 = std::sqrt(1 - rho_[i]) *
                (std::sqrt(1 - rho_[i] + N * rho_[i]) - std::sqrt(1 - rho_[i]));
    double c2 =
        std::pow(std::sqrt(1 - rho_[i] + N * rho_[i]) - std::sqrt(1 - rho_[i]),
                 2) *
        R1 / std::pow(N, 2);

    arma::mat mk = (1 - rho_[i]) * R + c1 * RJ2 + c2;

    lamk[i] = arma::eig_sym(mk);

    double tol = 1e-20;
    if (lamk[i].max() <= tol) {
      lamk[i] = arma::clamp(lamk[i], tol,
                            tol + std::numeric_limits<double>::epsilon());
    } else {
      lamk[i] = arma::clamp(lamk[i], tol, lamk[i].max());
    }

    pval[i] = SKAT_pval(Qw[i], lamk[i], tp.saddlepoint);
  }

  double pmin = pval.min();
  arma::vec qval = {0, 0, 0, 0, 0, 0, 0, 0};
  for (arma::uword i = 0; i < K - 1; i++) {
    qval[i] = Liu_qval_mod(pmin, lamk[i]);
  }

  arma::vec lam;
  bool success = arma::eig_sym(lam, R - (Rs * Rs.t()) / R1);
  if (!success) {
    return 1.;
  }
  lam = arma::abs(lam);

  if (lam.max() <= 0) {
    lam = arma::clamp(lam, 0, std::numeric_limits<double>::epsilon());
  } else {
    lam = arma::clamp(lam, 0, lam.max());
  }

  arma::vec tauk(K - 1);
  for (arma::uword i = 0; i < K - 1; i++) {
    tauk(i) = (1 - rho_[i]) * R2 / R1 + rho_[i] * R1;
  }
  double vp2 = 4. * (R3 / R1 - std::pow(R2, 2) / std::pow(R1, 2));
  double MuQ = arma::accu(lam);
  double VarQ = 2. * arma::accu(lam.t() * lam);
  double sd1 = std::sqrt(VarQ) / std::sqrt(VarQ + vp2);

  boost::math::chi_squared chisq(1);
  double q1 = boost::math::quantile(boost::math::complement(
      chisq, pmin > 0 ? pmin : std::sqrt(std::numeric_limits<double>::min())));
  double T0 = pmin;

  // Integration
  auto katint = [&](double xpar) -> double {
    double eta1 = std::numeric_limits<double>::max();
    for (arma::uword i = 0; i < K - 1; i++) {
      double val = (qval[i] - tauk[i] * xpar) / (1 - rho_[i]);
      if (val < eta1)
        eta1 = val;
    }
    double x = (eta1 - MuQ) * sd1 + MuQ;
    return SKAT_pval(x, lam, tp.saddlepoint) * boost::math::pdf(chisq, xpar);
  };

  double error_estimate;
  unsigned int max_depth = 5;
  double tolerance = 1e-9;
  double p_value = 1;
  // Can't calculate p-value, return alternate
  if (q1 < std::numeric_limits<double>::min() * 10) {
    return std::max(std::numeric_limits<double>::min(),
                    std::min(p_value, pmin * K));
  }
  p_value = T0 + boost::math::quadrature::gauss_kronrod<double, 15>::integrate(
                     katint, std::numeric_limits<double>::min() * 10, q1,
                     max_depth, tolerance, &error_estimate);

  if (p_value >= 1 || pmin >= 1) {
    std::cerr << "p_value: " << p_value << " pmin: " << pmin << "\n";
  }

  return std::max(std::numeric_limits<double>::min(),
                  std::min(p_value, pmin * K));
}

/**
 * @brief SKATC being SKATO but combining p-values with ACAT
 * @param gene A gene object representing the gene to be tested
 * @param transcript The transcript id to be tested
 * @param phenotypes The sample phenotypes to be tested
 * @param a The beta parameter for the beta distribution
 * @param b The beta parameter for the beta distribution
 * @param detail Whether to generate detailed output
 * @param linear Whether to fit a linear model
 * @return The p-value for the gene
 */
double Methods::SKATC(Gene &gene, arma::vec &phenotypes,
                      const std::string &transcript, int a, int b, bool detail,
                      bool linear) {
  if (linear) {
    lin_obj_->shuffle(phenotypes);
  } else {
    obj_->shuffle(phenotypes);
  }
  // ACAT p-value combination functions
  const double pi = boost::math::constants::pi<double>();
  auto acat = [&pi](arma::vec &p, arma::vec &w) -> double {
    return arma::accu(w % arma::tan((0.5 - p) * pi));
  };
  auto acat_p = [&pi](const double t, arma::vec &w) -> double {
    return 0.5 - atan(t / arma::accu(w)) / pi;
  };

#ifdef SKATCMC
  const double skatp =
      SKAT(gene, phenotypes, transcript, a, b, detail, linear, false);
  const double cmc1dfp = CMC1df(gene, phenotypes, transcript, false);

  arma::vec ps = {skatp, cmc1df};
  arma::vec ws = {1, 1};

  const double stat = acat(ps, ws);
  return acat_p(stat, ws);
#else
  const double skatp =
      SKAT(gene, phenotypes, transcript, a, b, detail, linear, false);
  const double burdenp = BURDEN(gene, phenotypes, transcript, false);

  arma::vec ps = {skatp, burdenp};
  arma::vec ws = {1, 1};

  const double stat = acat(ps, ws);
  return acat_p(stat, ws);
#endif
}

void Methods::check_weights(Gene &gene, const std::string &transcript, int a,
                            int b, bool no_weight) {
  if (no_weight) {
    arma::vec weights(gene.genotypes[transcript].n_cols, arma::fill::ones);
    gene.set_weights(transcript, weights);
    return;
  }
  const arma::sp_mat G(gene.genotypes[transcript]);
  arma::vec weights = gene.weights[transcript];

  if (kernel_ == Kernel::wLinear) {
    arma::vec maf(arma::mean(G, 0).t() / 2.);

    for (arma::uword i = 0; i < G.n_cols; i++) {
      weights(i) = std::pow(maf(i), a - 1) * std::pow(1 - maf(i), b - 1) /
                   boost::math::beta(a, b);
      // weights(i) = std::pow(maf(i), a - 1) * std::pow(1 - maf(i), b - 1);
    }
    gene.set_weights(transcript, weights);
  } else {
    gene.set_weights(transcript, weights);
  }
}

double Methods::call(Gene &gene, Covariates &cov, arma::vec &phenotypes,
                     const std::string &transcript, bool detail) {
  if (method_ == "BURDEN") {
    return BURDEN(gene, phenotypes, transcript, tp.nperm > 0);
  } else if (method_ == "CALPHA") {
    return CALPHA(gene, phenotypes, transcript);
  } else if (method_ == "CMC") {
    return CMC(gene, phenotypes, transcript, tp.cmcmaf);
  } else if (method_ == "CMC1df") {
    return CMC1df(gene, phenotypes, transcript, tp.nperm > 0);
  } else if (method_ == "RVT1") {
    return RVT1(gene, phenotypes, cov.get_covariate_matrix(), cov.get_coef(),
                transcript, tp.qtl);
  } else if (method_ == "RVT2") {
    return RVT2(gene, phenotypes, cov.get_covariate_matrix(), cov.get_coef(),
                transcript, tp.qtl);
  } else if (method_ == "SKAT") {
    return SKAT(gene, phenotypes, transcript, tp.a, tp.b, detail, tp.qtl,
                tp.nperm > 0);
  } else if (method_ == "SKATO") {
    return SKATO(gene, phenotypes, transcript, tp.a, tp.b, detail, tp.qtl);
  } else if (method_ == "SKATC") {
    return SKATC(gene, phenotypes, transcript, tp.a, tp.b, detail, tp.qtl);
  } else if (method_ == "VAAST") {
    return VAAST(gene, phenotypes, transcript, tp.vaast_site_penalty,
                 tp.group_size, detail, tp.biallelic, tp.soft_maf_filter,
                 tp.alternate_grouping);
  } else if (method_ == "VT") {
    return VT(gene, phenotypes, transcript);
  } else if (method_ == "WSS") {
    return WSS(gene, phenotypes, transcript);
  } else {
    throw(std::logic_error("Failed to find "));
  }
}
