//
// Created by Bohlender,Ryan James on 7/31/18.
//

#define ARMA_DONT_PRINT_ERRORS

#include <iomanip>
#include <cmath>
#include <chrono>

// Boost Math
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/roots.hpp>

#include <boost/format.hpp>

#include "methods.hpp"
#include "../data/gene.hpp"
#include "../data/covariates.hpp"
#include "vaast.hpp"

#include "../link/binomial.hpp"
#include "../link/gaussian.hpp"
#include "glm.hpp"
#include "fishertest.hpp"

constexpr std::array<double, 8> Methods::rho_;

// TODO Check with Chad / Yao re: replacing NAN values for MGIT
arma::vec rank(arma::vec &v, const char *direction) {
  if (strcmp(direction, "ascend") != 0 && strcmp(direction, "descend") != 0)
    throw (std::logic_error("Order argument for rank() must be either 'ascend' or 'descend'"));

  arma::uvec sort_indices;
  try {
    sort_indices = arma::sort_index(v, direction);
  } catch (const std::logic_error &e) {
    if (strcmp(direction, "ascend") == 0) {
      std::cerr << "NANs among ranked values. Replacing with 1.\n";
      v.replace(arma::datum::nan, 1);
    } else {
      std::cerr << "NANs among ranked values. Replacing with 0.\n";
      v.replace(arma::datum::nan, 0);
    }
    sort_indices = arma::sort_index(v, direction);
  }
  arma::vec sorted = v(sort_indices);

  arma::vec ranks = arma::vec(v.n_rows, arma::fill::zeros);
  arma::uword i = 0, j = 0;

  while (i < v.n_rows) {
    j = i + 1;
    // Find the next different value
    while (j < v.n_rows) {
      if (sorted(i) != sorted(j))
        break;
      j++;
    }
    // Adjusted rank
    for (arma::uword k = i; k <= j - 1; k++) {
      ranks(sort_indices(k)) = 1. + (i + j - 1.) / 2.0f;
    }
    // Update i
    i = j;
  }

  return ranks;
}

Methods::Methods(std::string method)
    : method_(std::move(method)),
      kernel_(Kernel::Linear) {
}

Methods::Methods(TaskParams &tp, const std::shared_ptr<Covariates> &cov)
    : method_(tp.method), tp_(tp) {
  if (tp.kernel == "Linear") {
    kernel_ = Kernel::Linear;
  } else if (tp.kernel == "wLinear") {
    kernel_ = Kernel::wLinear;
  }

  if ((tp.method == "SKAT" || tp.method == "SKATO" || tp.method == "BURDEN") && !tp.linear) {
    obj_ = std::make_shared<SKATR_Null>(cov);
    lin_obj_ = nullptr;
  } else if ((tp.method == "SKAT" || tp.method == "SKATO" || tp.method == "BURDEN") && tp.linear) {
    obj_ = nullptr;
    lin_obj_ = std::make_shared<SKATR_Linear_Null>(*cov);
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

double Methods::BURDEN(Gene &gene, const std::string &k, arma::vec &phenotypes, int a, int b) {
  obj_->shuffle(phenotypes);

  arma::sp_mat G(gene.get_matrix(k));

  check_weights(gene, k, a, b);

  arma::mat W = arma::diagmat(gene.get_weights(k));

  return std::pow(arma::accu(arma::sum(arma::diagmat(obj_->get_U0()) * G) * W), 2);
}

double Methods::CALPHA(Gene &gene, arma::vec &Y, const std::string &k) {
  arma::sp_mat X(gene.get_matrix(k));

  double nA = arma::sum(Y); // Case count
  double nU = Y.n_rows - nA; // Control count

  double p0 = nA / (nA + nU);

  arma::vec n(X.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i < X.n_cols; i++) {
    for (const auto &v: X.col(i)) {
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

double Methods::CMC(Gene &gene, arma::vec &Y, const std::string &k, double maf) {
  arma::mat X(gene.get_matrix(k));

  double N = Y.n_rows;
  double nA = arma::sum(Y);     // Case count
  double nU = N - nA;               // Control count

  arma::rowvec MAF = arma::mean(X, 0) / 2;

  // Collapse rare variants
  arma::uvec rare = arma::find(MAF < maf);
  arma::uvec common = arma::find(MAF >= maf);

  arma::mat Xnew;
  if (rare.size() <= 1) {
    Xnew = X;
  } else {
    arma::mat Xcollapse = arma::sum(X.cols(rare), 1);
    Xcollapse(arma::find(Xcollapse > 1)).ones();
    Xnew = X.cols(common);
    Xnew.insert_cols(Xnew.n_cols, Xcollapse);
  }

  if (tp_.nperm < 0) {
    // Rescale to -1, 0, 1
    Xnew -= 1;

    // Calculate two-sample Hotelling's T2 statistic
    arma::mat Xx = Xnew.rows(arma::find(Y == 1));
    arma::mat Yy = Xnew.rows(arma::find(Y == 0));

    arma::rowvec Xxmean = arma::mean(Xx);
    arma::rowvec Yymean = arma::mean(Yy);

    return arma::norm(Xxmean - Yymean);
  } else {
    // Rescale to -1, 0, 1
    Xnew -= 1;

    // Calculate two-sample Hotelling's T2 statistic
    arma::mat Xx = Xnew.rows(arma::find(Y == 1));
    arma::mat Yy = Xnew.rows(arma::find(Y == 0));

    arma::rowvec Xxmean = arma::mean(Xx);
    arma::rowvec Yymean = arma::mean(Yy);

    arma::mat COV = ((nA - 1.) * arma::cov(Xx) + (nU - 1.) * arma::cov(Yy)) / (N - 2.);
    arma::mat INV;
    if (!arma::inv_sympd(INV, COV)) {
      arma::pinv(INV, COV);
    }
    arma::mat ret = (Xxmean - Yymean) * INV * (Xxmean - Yymean).t() * nA * nU / N;
    auto p = static_cast<double>(Xxmean.n_elem);
    double stat = arma::as_scalar(ret) * (nA + nU - 1 - p) / (p * (nA + nU - 2)); // F(N, nA + nU - 1 - N) distributed
    if (stat < 0)
      stat = 0;
    if(tp_.nperm > 0) {
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
  }
}

double Methods::CMC1df(Gene &gene, arma::vec &Y, const std::string &k) {
  // Runtime for MDA OV with just fisher test and 10000 perms = 6544.95
  // Runtime for fast path with 10000 perms = 267.874
  if(tp_.nperm > 0) {
    arma::vec X(arma::sum(gene.get_matrix(k), 1));
    X(arma::find(X > 0)).ones();

    arma::uword ncase = 2 * arma::accu(Y);
    arma::uword ncont = 2 * arma::accu(1 - Y);

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

    return case_alt * cont_ref / (cont_alt * case_ref);
  } else {
    FisherTest fisherTest(gene, Y, k);
    return fisherTest.get_pval();
  }
}

double Methods::RVT1(Gene &gene, arma::vec &Y, arma::mat design, arma::vec &initial_beta, const std::string &k, bool linear) {
  // Runtime 100 perms naive initialization on macbook pro, ovarian data -- 2421.09
  // Runtime 100 perms prior initialization on macbook pro, ovarian data -- 1863.08 -- Poor initialization in permutation
  // Runtime 100 perms single prior initialization on macbook pro, ovarian data -- 1757.13
  if (linear) {
    // Quantitative trait
    arma::sp_mat X = gene.get_matrix(k).t();
    Gaussian link("identity");
    GLM<Gaussian> fit1(design, Y, link);
    arma::mat d2 = arma::join_horiz(design, arma::rowvec(arma::sum(X) / X.n_rows).t());
    GLM<Gaussian> fit2(d2, Y, link, fit1.beta_);

    double n = Y.n_elem;

    boost::math::chi_squared chisq(1);
    // TODO: Should be rank not n_rows
    double stat = (fit1.dev_ - fit2.dev_) / (fit2.dev_ / (n - d2.n_rows));
    if (stat < 0) {
      stat = std::numeric_limits<double>::epsilon();
    }
    return boost::math::cdf(boost::math::complement(chisq, stat));
  } else {
    // Binary trait
    arma::sp_mat X = gene.get_matrix(k).t();
    Binomial link("logit");
    GLM<Binomial> fit1(design, Y, link);
    arma::mat d2 = arma::join_horiz(design, arma::rowvec(arma::sum(X) / X.n_rows).t());
    GLM<Binomial> fit2(d2, Y, link, fit1.beta_);

    boost::math::chi_squared chisq(1);
    double stat = fit1.dev_ - fit2.dev_;
    if (stat < 0) {
      stat = std::numeric_limits<double>::epsilon();
    }
    return boost::math::cdf(boost::math::complement(chisq, stat));
  }
}

double Methods::RVT2(Gene &gene,
                     arma::vec &Y,
                     arma::mat design,
                     arma::vec &initial_beta,
                     const std::string &k,
                     bool linear) {
  if (linear) {
    // Quantitative trait
    arma::sp_mat X = gene.get_matrix(k).t();
    Gaussian link("identity");
    GLM<Gaussian> fit1(design, Y, link);
    arma::rowvec r = arma::conv_to<arma::rowvec>::from(arma::rowvec(arma::sum(X)) > 0);
    arma::mat d2 = arma::join_horiz(design, r.t());
    GLM<Gaussian> fit2(d2, Y, link, fit1.beta_);

    double n = Y.n_elem;

    boost::math::chi_squared chisq(1);
    double stat = (fit1.dev_ - fit2.dev_) / (fit2.dev_ / (n - d2.n_rows));
    if (stat < 0) {
      stat = std::numeric_limits<double>::epsilon();
    }
    return boost::math::cdf(boost::math::complement(chisq, stat));
  } else {
    // Binary trait
    arma::sp_mat X = gene.get_matrix(k).t();
    Binomial link("logit");
    GLM<Binomial> fit1(design, Y, link);
    arma::rowvec r = arma::conv_to<arma::rowvec>::from(arma::rowvec(arma::sum(X)) > 0);
    arma::mat d2 = arma::join_horiz(design, r.t());
    GLM<Binomial> fit2(d2, Y, link, fit1.beta_);

    boost::math::chi_squared chisq(1);
    double stat = fit1.dev_ - fit2.dev_;
    if (stat < 0) {
      stat = std::numeric_limits<double>::epsilon();
    }
    return boost::math::cdf(boost::math::complement(chisq, stat));
  }
}

double Methods::VAAST(Gene &gene,
					  arma::vec &Y,
					  const std::string &k,
					  double site_penalty,
					  arma::uword group_threshold,
					  bool detail,
					  bool biallelic,
					  double control_freq_cutoff,
					  bool legacy) {
  VAASTLogic vaast_logic
	  (gene,
	   Y,
	   k,
	   site_penalty,
	   group_threshold,
	   detail,
	   biallelic,
	   control_freq_cutoff,
	   legacy);
  return vaast_logic.score;
}

double Methods::VT(Gene &gene, const std::string &k, arma::vec &phenotypes) {
  // Convert data to match their format
  arma::vec pheno = arma::repmat(phenotypes, gene.get_matrix(k).n_cols, 1);
  if (!vt_obj_->is_initialized(k)) {
    vt_obj_->initialize(gene, pheno, k);
  }
  arma::vec phenoCount = pheno % vt_obj_->get_mCount(k); // Changes under permutation
  arma::vec csPhenoCount =
      arma::cumsum(vt_obj_->sum_groups(phenoCount, vt_obj_->get_oneToLen(k), k)); // Changes under permutation

  return arma::max((csPhenoCount - vt_obj_->get_csCountMeanpheno(k)) / vt_obj_->get_sqrtCsCountSquare(k));
}

double Methods::WSS(Gene &gene, arma::vec &Y, const std::string &k) {
  arma::mat X(gene.get_matrix(k));

  double nA = arma::sum(Y); // Case count
  double n = Y.n_rows;

  arma::vec mU = arma::sum(X, 0).t();
  arma::vec q = (mU + 1.) / (2. * n + 2.);

  arma::vec w = 1. / arma::sqrt(q % (1. - q));
  w.replace(arma::datum::nan, 0);

  arma::mat gamma_mat = X * arma::diagmat(w);
  gamma_mat.replace(arma::datum::nan, 0);

  arma::vec gamma = arma::sum(gamma_mat, 1);

  return arma::accu(gamma % Y);
}

std::string Methods::str() {
  return method_;
}

/**
 * @brief Calculate SKAT with p-value following Wu, Guan, Pankow (2017)
 * @param gene
 * @param cov
 * @param weights
 * @param k
 * @param phenotypes
 * @param a
 * @param b
 * @return
 */
double Methods::SKAT(Gene &gene,
                     const std::string &k,
                     arma::vec &phenotypes,
                     int a,
                     int b,
                     bool detail,
                     bool linear,
                     bool permute,
                     bool shuffle) {
  arma::sp_mat G(gene.get_matrix(k));

  if (shuffle) {
    if (linear) {
      lin_obj_->shuffle(phenotypes);
    } else {
      obj_->shuffle(phenotypes);
    }
  }

  check_weights(gene, k, a, b);
  arma::vec weights = gene.get_weights(k);

  arma::mat W = arma::diagmat(weights);
  // We're permuting, only calculate the Q-value
  if (permute) {
    arma::rowvec Zs;
    if (linear) {
      Zs = arma::sum(arma::diagmat(lin_obj_->get_U0()) * G) / std::sqrt(lin_obj_->get_s2());
    } else {
      Zs = arma::sum(arma::diagmat(obj_->get_U0()) * G);
    }
    arma::mat Z = Zs * W;

    double Q = arma::accu(arma::pow(Z, 2));

    return Q;
    // We're not permuting, return asymptotic p-values
  } else {

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
      Zs = arma::sum(arma::diagmat(lin_obj_->get_U0()) * G) / std::sqrt(lin_obj_->get_s2());
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
      gene.set_scores(k, variant_scores);
    }

    return SKAT_pval(Q, s);
  }
}

double Methods::SKATO(Gene &gene,
                      const std::string &k,
                      arma::vec &phenotypes,
                      int a,
                      int b,
                      bool detail,
                      bool linear) {
  if (linear) {
    lin_obj_->shuffle(phenotypes);
  } else {
    obj_->shuffle(phenotypes);
  }

  arma::sp_mat G(gene.get_matrix(k));
  arma::uword N = G.n_cols; // Variant count

  check_weights(gene, k, a, b);
  arma::vec weights = gene.get_weights(k);

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
    Zs = arma::sum(arma::diagmat(lin_obj_->get_U0()) * G) / std::sqrt(lin_obj_->get_s2());
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
    double c1 = std::sqrt(1 - rho_[i]) * (std::sqrt(1 - rho_[i] + N * rho_[i]) - std::sqrt(1 - rho_[i]));
    double c2 = std::pow(std::sqrt(1 - rho_[i] + N * rho_[i]) - std::sqrt(1 - rho_[i]), 2) * R1 / std::pow(N, 2);

    arma::mat mk = (1 - rho_[i]) * R + c1 * RJ2 + c2;

    lamk[i] = arma::eig_sym(mk);

    double tol = 1e-20;
    if (lamk[i].max() <= tol) {
      lamk[i] = arma::clamp(lamk[i], tol, tol + std::numeric_limits<double>::epsilon());
    } else {
      lamk[i] = arma::clamp(lamk[i], tol, lamk[i].max());
    }

    pval[i] = SKAT_pval(Qw[i], lamk[i]);
  }

  double pmin = pval.min();
  arma::vec qval = {0, 0, 0, 0, 0, 0, 0, 0};
  for (arma::uword i = 0; i < K - 1; i++) {
    qval[i] = Liu_qval_mod(pmin, lamk[i]);
  }

  arma::vec lam;
  arma::eig_sym(lam, R - (Rs * Rs.t()) / R1);
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
  double q1 = boost::math::quantile(boost::math::complement(chisq,
                                                            pmin > 0 ? pmin
                                                                     : std::sqrt(std::numeric_limits<double>::min())));
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
    return SKAT_pval(x, lam) * boost::math::pdf(chisq, xpar);
  };

  double error_estimate;
  unsigned int max_depth = 5;
  double tolerance = 1e-9;
  double p_value = 1;
  // Can't calculate p-value, return alternate
  if (q1 < std::numeric_limits<double>::min() * 10) {
    return std::max(std::numeric_limits<double>::min(), std::min(p_value, pmin * K));
  }
  p_value = T0 + boost::math::quadrature::gauss_kronrod<double, 15>::integrate(katint,
                                                                               std::numeric_limits<double>::min()
                                                                                   * 10,
                                                                               q1,
                                                                               max_depth,
                                                                               tolerance,
                                                                               &error_estimate);

  if (p_value >= 1 || pmin >= 1) {
    std::cerr << "p_value: " << p_value << " pmin: " << pmin << "\n";
  }

  return std::max(std::numeric_limits<double>::min(), std::min(p_value, pmin * K));
}

void Methods::check_weights(Gene &gene, const std::string &k, int a, int b, bool no_weight) {
  if (gene.is_weighted(k)) {
    return;
  } else if (no_weight) {
    arma::vec weights(gene.get_matrix(k).n_cols, arma::fill::ones);
    gene.set_weights(k, weights);
    return;
  }
  arma::mat G(gene.get_matrix(k));
  arma::vec weights(G.n_cols, arma::fill::ones);

  if (kernel_ == Kernel::wLinear) {
    arma::vec maf = arma::mean(G, 0).t() / 2.;

    for (arma::uword i = 0; i < G.n_cols; i++) {
      weights(i) = std::pow(maf(i), a - 1) * std::pow(1 - maf(i), b - 1) / boost::math::beta(a, b);
      //weights(i) = std::pow(maf(i), a - 1) * std::pow(1 - maf(i), b - 1);
    }
    if (method_ == "VAAST") {
      weights.replace(0, std::sqrt(std::numeric_limits<double>::min()));
    }
    gene.set_weights(k, weights);
  } else {
    gene.set_weights(k, weights);
  }
}

