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
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>

#include <boost/format.hpp>

#include "methods.hpp"
#include "../data/gene.hpp"
#include "../data/covariates.hpp"
#include "skat_adjust.hpp"
#include "../third_party/QFC/qfc2.hpp"
#include "vaast.hpp"

constexpr double Methods::rho_[];

// TODO Check with Chad / Yao re: replacing NAN values for MGIT
arma::vec rank(arma::vec &v, const char *direction) {
  if (strcmp(direction, "ascend") != 0 && strcmp(direction, "descend") != 0)
	throw (std::logic_error("Order argument must be either 'ascend' or 'descend'"));

  arma::uvec sort_indices;
  try {
	sort_indices = arma::sort_index(v, direction);
  } catch (const std::logic_error &e) {
	if (strcmp(direction, "ascend") == 0) {
	  std::cerr << "NANs among ranked values. Replacing with 1.\n";
	  v.replace(arma::datum::nan, 1);
	} else {
#ifndef NDEBUG
	  std::cerr << v.t();
	  std::cerr << e.what();
#endif
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

Methods::Methods(TaskParams &tp, Covariates &cov)
	: method_(tp.method) {
  if (tp.kernel == "Linear") {
	kernel_ = Kernel::Linear;
  } else if (tp.kernel == "wLinear") {
	kernel_ = Kernel::wLinear;
  } else if (tp.kernel == "IBS") {
	kernel_ = Kernel::IBS;
  } else if (tp.kernel == "wIBS") {
	kernel_ = Kernel::wIBS;
  } else if (tp.kernel == "Quadratic") {
	kernel_ = Kernel::Quadratic;
  } else if (tp.kernel == "twoWayX") {
	kernel_ = Kernel::twoWayX;
  }

  if (tp.alternate_permutation && !tp.linear) {
	obj_ = std::make_shared<SKATR_Null>(cov);
	lin_obj_ = nullptr;
  } else if(tp.alternate_permutation && tp.linear) {
    obj_ = nullptr;
    lin_obj_ = std::make_shared<SKATR_Linear_Null>(cov);
  }
}

/**
 * @brief Reset kernel to free memory.
 */
void Methods::clear(std::vector<std::string> &v) {
  for (const auto &k : v) {
	K_[k].reset();
	Q_sim_all[k].reset();
  }
  Q_sim_all.clear();
  re2.reset();
}

double Methods::BURDEN(Gene &gene, const std::string &k, bool shuffle, int a, int b) {
  if (shuffle) {
	obj_->shuffle();
  }

  arma::sp_mat G(gene.get_matrix(k));

  check_weights(gene, k, a, b);

  arma::mat W = arma::diagmat(gene.get_weights(k));

  return std::pow(arma::accu(arma::sum(arma::diagmat(obj_->get_U0()) * G) * W), 2);
}

double Methods::CALPHA(Gene &gene, Covariates &cov, const std::string &k) {
  arma::sp_mat X(gene.get_matrix(k));
  arma::vec Y(cov.get_phenotype_vector());

  double nA = arma::sum(Y); // Case count
  double nU = Y.n_rows - nA; // Control count

  double p0 = nA / (nA + nU);

  arma::vec n(X.n_cols, arma::fill::zeros);
  for(arma::uword i = 0; i < X.n_cols; i++) {
    for(const auto &v: X.col(i)) {
      if(v > 0)
        n(i)++;
    }
  }

  arma::vec g(X.n_cols, arma::fill::zeros);
  arma::uvec case_idx = arma::find(Y == 1);
  for(auto it = X.begin(); it != X.end(); it++) {
    if(arma::find(case_idx == it.row()).eval().n_elem > 0) {
      g(it.col())++;
    }
  }
  // arma::vec g = arma::sum(X.rows(arma::find(Y == 1)) > 0, 0).t();

  // Test statistic
  return arma::sum(arma::pow(g - (n * p0), 2) - (n * p0 * (1 - p0)));
}

double Methods::CMC(Gene &gene, Covariates &cov, const std::string &k, double maf) {
  arma::mat X(gene.get_matrix(k));
  arma::vec Y = cov.get_phenotype_vector();

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
	// Inversion failed
	arma::vec eigvals;
	arma::mat eigvecs;
	arma::eig_sym(eigvals, eigvecs, COV);

	eigvals = 1. / eigvals;
	arma::mat eigvecs_inv = arma::inv(eigvecs);

	INV = eigvecs * arma::diagmat(eigvals) * eigvecs_inv;
  }
  arma::mat ret = (Xxmean - Yymean) * INV * (Xxmean - Yymean).t() * nA * nU / N;
  return arma::as_scalar(ret);
}

#if 0
double Methods::SKAT(arma::mat &Xmat,
					 Covariates &cov,
					 arma::vec &weights,
					 const std::string &k,
					 bool shuffle,
					 int a,
					 int b) {
  arma::vec &Yvec = cov.get_phenotype_vector();

  // Randomize indices
  if (shuffle) {
	// arma::uvec ma_carriers = arma::find( arma::sum(Xmat) > 0 );
	cov.shuffle();
  }

  // Check for weighted kernel
  if (kernel_ == Kernel::wIBS || kernel_ == Kernel::wLinear) {
	// Check for weights if kernel hasn't been generated
	if (K_[k].n_rows == 0) {
	  weights.reshape(Xmat.n_cols, 1);
	  arma::vec maf = arma::mean(Xmat, 0).t() / 2.;

	  for (arma::uword i = 0; i < Xmat.n_cols; i++) {
		weights(i) = std::pow(maf(i), a - 1) * std::pow(1 - maf(i), b - 1) / boost::math::beta(a, b);
#if 0
		std::cerr << "a: " << a << " b: " << b << " ";
		std::cerr << "maf(i) = " << maf(i) << " weights(i) = " << weights(i) << "\n";
#endif
	  }
	}
  } else {
	if (K_[k].n_rows == 0) {
	  weights.ones(Xmat.n_cols, 1);
	}
  }

  arma::mat tXmat;
  arma::uword n = Xmat.n_cols; // Number of samples
  arma::uword p = Xmat.n_rows; // Number of variants

  // Only build kernel once
  if (K_[k].n_rows != n) {
	switch (kernel_) {
	case Kernel::Linear:K_[k] = kernel_Linear(Xmat);
	  break;
	case Kernel::wLinear:K_[k] = kernel_wLinear(Xmat, weights);
	  break;
	case Kernel::IBS:tXmat = Xmat.t();
	  K_[k] = kernel_IBS(tXmat, n, p);
	  break;
	case Kernel::wIBS:tXmat = Xmat.t();
	  K_[k] = kernel_wIBS(tXmat, n, p, weights);
	  break;
	case Kernel::Quadratic:K_[k] = kernel_Quadratic(Xmat);
	  break;
	case Kernel::twoWayX:tXmat = Xmat.t();
	  K_[k] = kernel_twoWayX(tXmat, n, p);
	  break;
	default:throw (std::runtime_error("Incorrect Kernel specification"));
	}
  }

  arma::vec Z = Yvec(cov.get_indices()) - cov.get_probability()(cov.get_indices());
  arma::mat ret = Z.t() * K_[k] * Z / 2;
  return ret(0, 0);
}

double Methods::SKATO(Gene &gene,
					  Covariates &cov,
					  arma::vec &weights,
					  const std::string &k,
					  bool shuffle,
					  int a,
					  int b,
					  bool adjust) {
  // Randomize indices
  if (shuffle) {
	// arma::uvec ma_carriers = arma::find( arma::sum(Z) > 0 );
	cov.shuffle();
  }

  arma::mat &Z = gene.get_matrix(k);

  if (Z.n_rows < 2000 && adjust) {
	// SKAT_Adjust(Gene &gene, Covariates &cov, const std::string &k, const std::string &kernel, int a, int b);
	std::string kernel = "Linear";
	if (kernel_ == Kernel::wLinear) {
	  kernel = "wLinear";
	}

	SKAT_Adjust skat_adjust(gene, cov, k, kernel, a, b, re2, Q_sim_all);

	return skat_adjust.p_value;
  }
  // Fallthrough, no small sample size adjustment

  arma::uword pm = Z.n_cols;
  arma::uword nr = 8;

  arma::vec rall(nr);

  for (arma::uword i = 0; i < nr; i++) {
	if (rho_[i] > 0.999) {
	  rall(i) = 0.999;
	} else {
	  rall(i) = rho_[i];
	}
  }

  // SKAT_Optimal_Logistic
  arma::vec res = cov.get_phenotype_vector()(cov.get_indices()) - cov.get_probability()(cov.get_indices());
  arma::mat X1 = cov.get_covariate_matrix().t();

  arma::vec pi_1 = cov.get_probability()(cov.get_indices()) % (1 - cov.get_probability()(cov.get_indices()));
#if 1
  arma::mat Z1 = (arma::diagmat(arma::sqrt(pi_1)) * Z
	  - (arma::diagmat(arma::sqrt(pi_1)) * X1) * arma::inv_sympd(X1.t() * (arma::diagmat(pi_1) * X1))
		  * (X1.t() * (arma::diagmat(pi_1) * Z)))
	  / arma::datum::sqrt2;
#else
  // Schur product is slow -- Roughly 25% speed improvement using diagonal matrices
  arma::mat Z1 = ((Z.each_col() % arma::sqrt(pi_1))
	  - (X1.each_col() % arma::sqrt(pi_1)) * arma::inv_sympd(X1.t() * (X1.each_col() % pi_1))
		  * (X1.t() * (Z.each_col() % pi_1)))
	  / arma::datum::sqrt2;
#endif
  arma::mat temp_res_out;
  SKAT_Optimal_GetQ skat_optimal_getQ(Z, res, rall, temp_res_out, 0);

  std::vector<arma::vec> lambda_all;
  std::vector<LiuParam> liu_param_all;

  for (arma::uword i = 0; i < nr; i++) {
	double rcorr = rho_[i];
	if (rcorr == 1)
	  rcorr = 0.999; // Prevents cholesky decomposition failure
	arma::mat Rm = arma::diagmat(arma::vec(pm).fill(1 - rcorr)) + arma::mat(pm, pm).fill(rcorr);
	arma::mat Z2 = Z1 * arma::chol(Rm).t();
	arma::mat K1 = Z2.t() * Z2;

	arma::vec eigval;
	arma::mat eigvec;

	arma::eig_sym(eigval, eigvec, K1);

	arma::uvec IDX1 = arma::find(eigval >= 0); // Positive eigenvalues
	arma::uvec IDX2 = arma::find(eigval > arma::mean(eigval(IDX1)) / 100000);

	// eigval(IDX2) == their lambda.temp for a given i
	arma::vec lambda = eigval(IDX2);

	lambda_all.emplace_back(lambda);
  }

  // SKAT_Optimal_Each_Q
  SKAT_EachQ skatEachQ(skat_optimal_getQ.Q_r, lambda_all);
  SKATParam skp(Z1);

  double pval = 0;

  pval = SKAT_Optimal_Pvalue_Davies(skatEachQ.pmin_q, skp, rall, skatEachQ.pmin);
#ifndef NDEBUG
  std::cerr << "pval: " << pval << "\n";
#endif

  arma::vec pval_each;
  double multi = 3;

  pval_each = skatEachQ.pval;
  arma::uvec IDX = arma::find(pval_each > 0);

  double pval1 = arma::min(pval_each) * multi;
  if (pval <= 0 || IDX.n_rows < rall.n_rows) {
	pval = pval1;
  }

  if (pval == 0) {
	if (IDX.n_rows > 0) {
	  pval = arma::min(pval_each(IDX));
	}
  }

#ifndef NDEBUG
  std::cerr << "pval_each: " << pval_each.t();
  std::cerr << "qval_each: " << skat_optimal_getQ.Q_r.t();
#endif

  return pval;
}
#endif

double Methods::Vaast(Gene &gene,
					  Covariates &cov,
					  const std::string &k,
					  bool score_only_minor,
					  bool score_only_alternative,
					  double site_penalty,
					  arma::uword group_threshold,
					  bool detail) {
#if 1
  VAAST vaast(gene, cov, k, score_only_minor, score_only_alternative, site_penalty, group_threshold, detail);
  return vaast.score;
#else
  arma::mat Xmat = gene.get_matrix(k);
  arma::vec Yvec = cov.get_phenotype_vector();

  check_weights(gene, k, 1, 25, true);
  // Do Grouping
  // Grouping logic: Given a group size, collapse variants with adjacent CASM scores. Collapse LGD variants separately.
  arma::vec weights = gene.get_weights(k);

  arma::mat Xnew;
  arma::vec Wnew;
  std::vector<arma::uvec> groups;
  arma::uvec ungrouped;

  std::tie(Xnew, Wnew, groups, ungrouped) = variant_grouping(Xmat, Yvec, weights, group_threshold);

  double n_case = arma::sum(Yvec);
  double n_control = Yvec.n_rows - n_case;

  arma::vec case_allele1(Xnew.n_cols);
  arma::vec control_allele1(Xnew.n_cols);

  case_allele1 = Xnew.t() * Yvec;
  control_allele1 = Xnew.t() * (1. - Yvec);

  arma::vec case_allele0 = 2 * n_case - case_allele1;
  arma::vec control_allele0 = 2 * n_control - control_allele1;

  // Get ln likelihood  of each variant
  arma::vec log_lh = LRT(case_allele1, control_allele1, case_allele0, control_allele0);

  arma::vec vaast_site_scores = 2.0 * (log_lh + arma::log(Wnew)) - site_penalty;

  // Mask sites
  arma::uvec mask = variant_bitmask(vaast_site_scores,
									case_allele1,
									control_allele1,
									case_allele0,
									control_allele0,
									score_only_minor,
									score_only_alternative);

  vaast_site_scores(mask).zeros();

  double total_score = arma::accu(vaast_site_scores);

  // Expand the scores
  arma::vec expanded_scores;
  if (group_threshold == 1) {
	expanded_scores = vaast_site_scores;
  } else {
	expanded_scores.reshape(Xmat.n_cols, 1);
	arma::uword j = 0;
	for (const auto &v : groups) {
	  for (const auto &i : v) {
		// Recalculate group score without current variant
		arma::uvec index = arma::regspace<arma::uvec>(0, Xmat.n_cols - 1);
		index.shed_row(i);

		arma::mat altX;
		arma::vec altW;
		std::vector<arma::uvec> altgroups;
		arma::uvec altungrouped;

		std::tie(altX, altW, altgroups, altungrouped) =
			variant_grouping(Xmat.cols(index), Yvec, weights(index), group_threshold);

		arma::vec case_count1 = altX.t() * Yvec;
		arma::vec cont_count1 = altX.t() * (1. - Yvec);
		arma::vec case_count0 = 2 * n_case - case_count1;
		arma::vec cont_count0 = 2 * n_control - cont_count1;

		arma::uvec altmask = variant_bitmask(vaast_site_scores,
											 case_count1,
											 cont_count1,
											 case_count0,
											 cont_count0,
											 score_only_minor,
											 score_only_alternative);

		arma::vec new_group_loglh = LRT(case_count1, cont_count1, case_count0, cont_count0);
		arma::vec new_score = 2.0 * (new_group_loglh + arma::log(altW)) - site_penalty;

		new_score(altmask).zeros();

		expanded_scores(i) = (total_score - arma::accu(new_score) < 0) ? 0 : total_score - arma::accu(new_score);
	  }
	  j++;
	}
	expanded_scores(ungrouped) = vaast_site_scores(arma::span(j, j + ungrouped.n_rows - 1));
  }
  expanded_scores(arma::find(expanded_scores < 0)).zeros();

  // Store scores for detailed output
  if (detail) {
	gene.set_scores(k, expanded_scores);
  }

  return total_score;
#endif
}

double Methods::VT(Gene &gene, Covariates &cov, const std::string &k) {
  arma::mat X(gene.get_matrix(k));
  arma::vec Y = cov.get_phenotype_vector();

  // All variants should be the minor allele
  arma::vec maf = ((1. + arma::sum(X, 0)) / (2. + 2. * X.n_rows)).t();
  arma::vec hmaf = arma::unique(maf);

  arma::vec zscores(hmaf.n_rows - 1, arma::fill::zeros);
  arma::vec res = Y - cov.get_fitted();

  for (arma::uword i = 0; i < hmaf.n_rows - 1; i++) {
	arma::mat Xmat_subset = X.cols(arma::find(maf < hmaf[i + 1]));
	double znum = arma::sum(arma::sum(arma::diagmat(res) * Xmat_subset, 0));
	double zden = std::sqrt(arma::sum(arma::sum(arma::pow(X.cols(arma::find(maf < hmaf[i + 1])), 2))));
	zscores(i) = znum / zden;
  }

  try {
	return arma::max(zscores);
  } catch (std::logic_error &e) {
	// TODO Check for better failure condition in the original paper. What should the value be in the case of a single variant?
	return 0;
  }
}

double Methods::WSS(Gene &gene, Covariates &cov, const std::string &k) {
  arma::mat X(gene.get_matrix(k));
  arma::vec Y = cov.get_phenotype_vector();

  double nA = arma::sum(Y); // Case count
  double nU = Y.n_rows - nA; // Control count
  double n = Y.n_rows;

  arma::vec mU = arma::sum(X.rows(arma::find(Y == 0)), 0).t();
  arma::vec q = (mU + 1.) / (2. * nU + 2.);

  arma::mat w = arma::diagmat(1. / arma::sqrt(n * (arma::diagmat(q) * (1. - q))));

  arma::mat gamma_mat = X * w;
  gamma_mat.replace(arma::datum::nan, 0);

  arma::vec gamma = arma::sum(gamma_mat, 1);

  arma::vec ranks = rank(gamma, "ascend");

  return arma::sum(ranks(arma::find(Y > 0)));
}

std::string Methods::str() {
  return method_;
}

/*
 * Kernel Member Functions
 */

arma::mat Methods::kernel_Linear(arma::mat &Xmat) {
  return Xmat * Xmat.t();
}

arma::mat Methods::kernel_wLinear(arma::mat &Xmat, arma::vec &weights) {
  return Xmat * arma::diagmat(weights) * arma::diagmat(weights) * Xmat.t();
}

/**
 * @brief Identity-by-Share Kernel
 * @param Xmat The genotype matrix.
 * @param n The number of samples.
 * @param p The number of variants.
 * @return Kernel as matrix.
 */
arma::mat Methods::kernel_IBS(arma::mat &Xmat, arma::uword &n, arma::uword &p) {
  arma::mat K = arma::eye(n, n);

  // For kernel calculation
  double temp, diff;

  for (int i = 0; i < n - 1; i++) {
	for (int j = i + 1; j < n; j++) {
	  temp = 0;
	  for (int k = 0; k < p; k++) {
		diff = 2 - abs(Xmat(i * p + k) - Xmat(j * p + k));
		temp += diff;
	  }
	  K(i * n + j) = K(j * n + i) = temp / 2. / p;
	}
  }
  return K;
}

/**
 * @brief Weighted Identity-by-Share Kernel
 * @param Xmat Genotype matrix.
 * @param n The number of samples.
 * @param p The number of variants.
 * @param weights The weights for each variant.
 * @return Kernel as matrix.
 */
arma::mat Methods::kernel_wIBS(arma::mat &Xmat, arma::uword &n, arma::uword &p, arma::vec &weights) {
  double diff, temp;

  arma::mat K = arma::eye(n, n);
  double w_total = arma::sum(weights);

  for (int i = 0; i < n - 1; i++) {
	for (int j = i + 1; j < n; j++) {
	  temp = 0;
	  for (int k = 0; k < p; k++) {
		diff = std::abs(Xmat(i * p + k) - Xmat(j * p + k));
		temp += diff * weights(k);
	  }
	  K(i * n + j) = K(j * n + i) = 1 - temp / 2 / w_total;
	}
  }
  return K;
}

/**
 * @brief Quadratic Kernel
 * @param Xmat The genotype matrix.
 * @return Quadratic kernel as matrix.
 */
arma::mat Methods::kernel_Quadratic(arma::mat &Xmat) {
  return arma::pow(Xmat * Xmat.t() + 1, 2);
}

/**
 * @brief Two Way kernel
 * @param Xmat The genotype matrix
 * @param n The number of samples.
 * @param p The number of variants.
 * @return Two way kernel as matrix.
 */
arma::mat Methods::kernel_twoWayX(arma::mat &Xmat, arma::uword n, arma::uword p) {
  arma::mat K = arma::zeros(n, n);

  double temp, temp1, temp2;

  for (int i = 0; i < n; i++) {
	for (int j = i; j < n; j++) {
	  temp = 1;

	  for (int k = 0; k < p; k++) {
		temp2 = Xmat(i * p + k) * Xmat(j * p + k);
		temp += temp2;
		if (k == 0) {
		  temp1 = temp2;
		  continue;
		}

		temp += temp1 * Xmat(i * p + k) * Xmat(j * p + k);
		temp1 += temp2;

	  }
	  K(i * n + j) = K(j * n + i) = temp;
	}
  }
  return K;
}

/**
 * @brief Calculate SKAT with p-value following Wu, Guan, Pankow (2017)
 * @param gene
 * @param cov
 * @param weights
 * @param k
 * @param shuffle
 * @param a
 * @param b
 * @return
 */
double Methods::SKATR(Gene &gene, const std::string &k, bool shuffle, int a, int b, bool detail, bool linear) {
  arma::sp_mat G(gene.get_matrix(k));

  if (shuffle) {
	if (linear) {
	  lin_obj_->shuffle();
	} else {
	  obj_->shuffle();
	}
  }

  check_weights(gene, k);
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
  if(linear) {
    Gs = G.t() * G - tmp.t() * tmp;
    Zs = arma::sum(arma::diagmat(lin_obj_->get_U0()) * G) / std::sqrt(lin_obj_->get_s2());
  } else {
    Gs = (arma::diagmat(obj_->get_Yv()) * G).t() * G - tmp.t() * tmp;
    Zs = arma::sum(arma::diagmat(obj_->get_U0()) * G);
  }

  arma::mat R = (Gs * W) * W;
  arma::mat Z = Zs * W;

  arma::vec s;
  arma::mat U, V;
  arma::svd(U, s, V, R);

  double Q = arma::accu(arma::pow(Z, 2));

  if (detail) {
	arma::vec variant_scores = arma::sum(arma::pow(Z, 2), 0).t();
	gene.set_scores(k, variant_scores);
  }

  return SKAT_pval(Q, s);
}

double Methods::SKATRO(Gene &gene, const std::string &k, bool shuffle, int a, int b, bool detail, bool linear) {
  if (shuffle) {
	if (linear) {
	  lin_obj_->shuffle();
	} else {
	  obj_->shuffle();
	}
  }

  arma::mat G(gene.get_matrix(k));
  arma::uword N = G.n_cols; // Variant count

  check_weights(gene, k);
  arma::vec weights = gene.get_weights(k);

  arma::mat W = arma::diagmat(weights);

  arma::mat tmp;
  if(linear) {
    tmp = lin_obj_->get_Ux().t() * G;
  } else {
	tmp = obj_->get_Ux().t() * G;
  }


  arma::mat Gs;
  arma::rowvec Zs;
  if(linear) {
	Gs = G.t() * G - tmp.t() * tmp;
	Zs = arma::sum(arma::diagmat(lin_obj_->get_U0()) * G) / std::sqrt(lin_obj_->get_s2());
  } else {
	Gs = (arma::diagmat(obj_->get_Yv()) * G).t() * G - tmp.t() * tmp;
	Zs = arma::sum(arma::diagmat(obj_->get_U0()) * G);
  }

  arma::mat R = (Gs * W).t() * W;
  arma::mat Z = Zs * W;

  arma::vec s;
  arma::mat U, V;
  arma::svd(U, s, V, R);

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
  double R3 = arma::accu(Rs % arma::sum(R.each_col() % Rs).t());

  arma::mat RJ2(Rs.n_rows, Rs.n_rows, arma::fill::zeros);
  for (arma::uword i = 0; i < Rs.n_rows; i++) {
	RJ2.row(i) = (Rs(i) + Rs.t()) / N;
  } // Replacement for R's outer(Rs, Rs, '+')

  std::vector<arma::vec> lamk(K - 1);
  for (arma::uword i = 0; i < K; i++) {
	// Pure burden
	if (rho_[i] == 1) {
	  boost::math::chi_squared chisq(1); // 1-df chisq
	  pval[i] = 1 - boost::math::cdf(chisq, Qb / R1);
	  break;
	}

	// Setup Davies
	double c1 = std::sqrt(1 - rho_[i]) * (std::sqrt(1 - rho_[i] + N * rho_[i]) - std::sqrt(1 - rho_[i]));
	double c2 = std::pow(std::sqrt(1 - rho_[i] + N * rho_[i]) - std::sqrt(1 - rho_[i]), 2) * R1 / std::pow(N, 2);

	arma::mat mk = (1 - rho_[i]) * R + c1 * RJ2 + c2;

	arma::mat Utmp, Vtmp;
	arma::svd_econ(Utmp, lamk[i], Vtmp, mk);

	double tol = 1e-20;
	lamk[i] = arma::clamp(lamk[i], tol, lamk[i].max());

	pval[i] = SKAT_pval(Qw[i], lamk[i]);
  }

  double pmin = pval.min();
  arma::vec qval = {0, 0, 0, 0, 0, 0, 0, 0};
  for (arma::uword i = 0; i < K - 1; i++) {
	qval[i] = Liu_qval_mod(pmin, lamk[i]);
  }

  arma::vec lam;
  arma::mat Utmp, Vtmp;
  arma::svd(Utmp, lam, Vtmp, R - (Rs * Rs.t()) / R1);

  lam = arma::clamp(lam, 0, lam.max()).eval()(arma::span(0, lam.n_rows - 1), arma::span::all);

  arma::vec tauk(K - 1);
  for (arma::uword i = 0; i < K - 1; i++) {
	tauk(i) = (1 - rho_[i]) * R2 / R1 + rho_[i] * R1;
  }
  double vp2 = 4. * (R3 / R1 - (R2 * R2) / (R1 * R1));
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
  double tolerance = 1e-25;
  double p_value = T0 + boost::math::quadrature::gauss_kronrod<double, 21>::integrate(katint,
																					  std::numeric_limits<double>::min()
																						  * 10,
																					  q1,
																					  max_depth,
																					  tolerance,
																					  &error_estimate);

  if (p_value >= 1 || pmin >= 1) {
	std::cerr << "p_value: " << p_value << " pmin: " << pmin << "\n";
  }

  return std::min(p_value, pmin * K);
}

double Methods::Liu_qval_mod(double pval, arma::vec lambda) {
  arma::vec c1{
	  arma::accu(lambda),
	  arma::accu(arma::pow(lambda, 2)),
	  arma::accu(arma::pow(lambda, 3)),
	  arma::accu(arma::pow(lambda, 4))
  };

  double muQ = c1[0];
  double sigmaQ = std::sqrt(2 * c1[1]);

  double s1 = c1[2] / std::pow(c1[1], 3. / 2.);
  double s2 = c1[3] / std::pow(c1[1], 2);

  double beta1 = 2 * arma::datum::sqrt2 * s1;
  double beta2 = 12 * s2;

  double type1 = 0;
  double a, d, l;
  if (s1 * s1 > s2) {
	a = 1. / (s1 - std::sqrt(s1 * s1 - s2));
	d = s1 * std::pow(a, 3) - a * a;
	l = a * a - 2 * d;
  } else {
	type1 = 1;
	l = 1. / s2;
	a = std::sqrt(l);
	d = 0;
  }

  double muX = l * d;
  double sigmaX = arma::datum::sqrt2 * a;
  double df = l;

  boost::math::chi_squared chisq(df);
  double q = boost::math::quantile(boost::math::complement(chisq,
														   pval > 0 ? pval
																	: std::sqrt(std::numeric_limits<double>::min())));
  return (q - df) / std::sqrt(2 * df) * sigmaQ + muQ;
}

double Methods::Saddlepoint(double Q, arma::vec lambda) {
  // Check for valid input
  if (Q <= 0) {
	return 1;
  }

  double d = lambda.max();
  if (d == 0) {
	return Liu_pval(Q, lambda);
  }
  arma::vec ulambda = lambda / d;
  Q /= d;

  if (ulambda.has_nan()) {
	ulambda.replace(arma::datum::nan, 0);
	std::cerr << ulambda.t();
  }

  auto k0 = [&](double &zeta) -> double {
	return -arma::accu(arma::log(1 - 2 * (zeta * ulambda))) / 2;
  };
  auto kprime0 = [&](double &zeta) -> double {
	return arma::accu(ulambda / (1 - 2 * zeta * ulambda));
  };
  auto kpprime0 = [&](double &zeta) -> double {
	return 2 * arma::accu(arma::pow(ulambda, 2) / arma::pow(1 - 2 * (zeta * ulambda), 2));
  };
  auto hatzetafn = [&](double zeta) -> double {
	return kprime0(zeta) - Q;
  };

  arma::uword n = ulambda.size();

  double lmin, lmax;
  if (arma::any(ulambda < 0)) {
	lmin = arma::max(1 / (2 * ulambda(arma::find(ulambda < 0)))) * 0.99999;
  } else if (Q > sum(ulambda)) {
	lmin = -0.01;
  } else {
	lmin = -static_cast<int>(ulambda.size()) / (2. * Q);
  }
  lmax = arma::min(1 / (2 * ulambda(arma::find(ulambda > 0)))) * 0.99999;

  // Root finding
#if 0
  int digits = std::numeric_limits<double>::digits - 3;
  boost::math::tools::eps_tolerance<double> tol(digits);
  boost::uintmax_t max_iter = 1000;
  std::pair<double, double>
	  tmp = boost::math::tools::toms748_solve(hatzetafn, lmin, lmax, tol, max_iter);
#elif 0
  int digits = std::numeric_limits<double>::digits - 3;
  boost::math::tools::eps_tolerance<double> tol(digits);
  boost::uintmax_t max_iter = 1000;
  double factor = (lmax - lmin) / 100;
  std::pair<double, double>
	  tmp = boost::math::tools::bracket_and_solve_root(hatzetafn, lmin, factor, true, tol, max_iter);
#else
  int digits = std::numeric_limits<double>::digits - 3;
  boost::math::tools::eps_tolerance<double> tol(digits);
  boost::uintmax_t max_iter = 1000;
  std::pair<double, double>
	  tmp = boost::math::tools::bisect(hatzetafn, lmin, lmax, tol, max_iter);

#endif

  double hatzeta = tmp.first + (tmp.second - tmp.first) / 2.;

  double w = sgn(hatzeta) * std::sqrt(2 * (hatzeta * Q - k0(hatzeta)));
  double v = hatzeta * std::sqrt(kpprime0(hatzeta));

  if (std::abs(hatzeta) < 1e-4 || std::isnan(w) || std::isnan(v)) {
	return Liu_pval(Q, lambda);
  } else {
	boost::math::normal norm(0, 1);
	return boost::math::cdf(boost::math::complement(norm, w + log(v / w) / w));
  }
}

template<class T>
int Methods::sgn(T x) {
  return (T(0) < x) - (x < T(0));
}

double Methods::Liu_pval(double Q, arma::vec lambda) {
  arma::vec c1{
	  arma::accu(lambda),
	  arma::accu(arma::pow(lambda, 2)),
	  arma::accu(arma::pow(lambda, 3)),
	  arma::accu(arma::pow(lambda, 4))
  };
  double muQ = c1[0];
  double sigmaQ = std::sqrt(2 * c1[1]);

  double s1 = c1[2] / std::pow(c1[1], 1.5);
  double s2 = c1[3] / std::pow(c1[1], 2);

  double a, d, l;
  if (s1 * s1 > s2) {
	a = 1. / (s1 - std::sqrt(s1 * s1 - s2));
	d = s1 * std::pow(a, 3) - a * a;
	l = a * a - 2 * d;
  } else {
	l = 1. / s2;
	a = std::sqrt(l);
	d = 0;
  }
  double muX = l * d;
  double sigmaX = arma::datum::sqrt2 * a;
  double df = l;

  double Qnorm = (Q - muQ) / sigmaQ * sigmaX + muX;

  try {
	boost::math::non_central_chi_squared chisq(df, d);
	return 1 - boost::math::cdf(chisq, Qnorm);
  } catch (std::exception &e) {
	return 1;
  }
}

double Methods::SKAT_pval(double Q, arma::vec lambda) {
  if (std::isnan(Q)) {
	return 1;
  }
  std::vector<double> lb1 = arma::conv_to<std::vector<double>>::from(lambda);
  std::vector<double> nc1(lb1.size(), 0);
  std::vector<int> df(lb1.size(), 1);
  double sigma = 0;
  int lim1 = 1000000;
  double acc = 1e-9;

  QFC qfc(lb1, nc1, df, sigma, Q, lim1, acc);

  double pval = 1 - qfc.get_res();
  if (pval >= 1 || pval <= 0 || qfc.get_fault() > 0) {
	pval = Saddlepoint(Q, lambda);
  }
  return pval;
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
	}
	if (method_ == "VAAST") {
	  weights.replace(0, std::sqrt(std::numeric_limits<double>::min()));
	}
	gene.set_weights(k, weights);
  } else {
	gene.set_weights(k, weights);
  }
}

