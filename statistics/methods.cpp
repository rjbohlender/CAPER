//
// Created by Bohlender,Ryan James on 7/31/18.
//

#define ARMA_DONT_PRINT_ERRORS

#include <iomanip>
#include <cmath>
// Boost Math
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/chi_squared.hpp>

// Boost Integration
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>

#include "methods.hpp"
#include "../data/gene.hpp"
#include "../data/covariates.hpp"

extern "C" {
#include "../third_party/SKAT/skat.hpp"
};

#include "../third_party/SKAT/qfc2.hpp"

constexpr double SKATParam::rho_[];
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

Methods::Methods(std::string method, std::string kernel)
	: method_(std::move(method)) {
  if (kernel == "Linear") {
	kernel_ = Kernel::Linear;
  } else if (kernel == "wLinear") {
	kernel_ = Kernel::wLinear;
  } else if (kernel == "IBS") {
	kernel_ = Kernel::IBS;
  } else if (kernel == "wIBS") {
	kernel_ = Kernel::wIBS;
  } else if (kernel == "Quadratic") {
	kernel_ = Kernel::Quadratic;
  } else if (kernel == "twoWayX") {
	kernel_ = Kernel::twoWayX;
  }
}

/**
 * @brief Reset kernel to free memory.
 */
void Methods::clear(std::vector<std::string> &v) {
  for (const auto &k : v) {
	K_[k].reset();
  }
}

double Methods::BURDEN(arma::mat &Xmat, Covariates &cov, arma::vec &weights) {
  arma::vec &Yvec = cov.get_phenotype_vector();
  arma::vec res = Yvec(cov.get_indices()) - cov.get_probability()(cov.get_indices());

  arma::rowvec S = res.t() * Xmat;

  return std::pow(arma::sum(S), 2);
}

double Methods::CALPHA(arma::mat &Xmat, arma::vec &Yvec) {
  double nA = arma::sum(Yvec); // Case count
  double nU = Yvec.n_rows - nA; // Control count

  double p0 = nA / (nA + nU);

  arma::vec n = arma::conv_to<arma::colvec>::from(arma::sum(Xmat > 0, 0));
  arma::vec g = arma::conv_to<arma::colvec>::from(arma::sum(Xmat.rows(arma::find(Yvec == 1)) > 0, 0));

  // Test statistic
  return arma::sum(arma::pow(g - (n * p0), 2) - (n * p0 * (1 - p0)));
}

double Methods::CMC(arma::mat &Xmat, arma::vec &Yvec, double maf) {
  int N = static_cast<int>(Yvec.n_rows);
  int nA = static_cast<int>(arma::sum(Yvec));
  int nU = N - nA;

  arma::rowvec MAF = arma::mean(Xmat, 0) / 2;

  // Collapse rare variants
  arma::uvec rare = arma::find(MAF < maf);
  arma::uvec common = arma::find(MAF >= maf);

  arma::mat Xnew;
  if (rare.size() <= 1) {
	Xnew = Xmat;
  } else {
	arma::mat Xcollapse = arma::sum(Xmat.cols(rare), 1);
	Xcollapse(arma::find(Xcollapse > 1)).ones();
	Xnew = Xmat.cols(common);
	Xnew.insert_cols(Xnew.n_cols, Xcollapse);
  }

  // Rescale to -1, 0, 1
  Xnew -= 1;

  arma::mat Xx = Xnew.rows(arma::find(Yvec == 1));
  arma::mat Yy = Xnew.rows(arma::find(Yvec == 0));

  arma::rowvec Xxmean = arma::mean(Xx, 0);
  arma::rowvec Yymean = arma::mean(Yy, 0);

  // Center matrices
  arma::mat Dx = Xx;
  arma::mat Dy = Yy;

  for (arma::uword i = 0; i < Dx.n_cols; i++) {
	Dx.col(i) -= Xxmean(i);
	Dy.col(i) -= Yymean(i);
  }

  arma::mat COV = (Dx.t() * Dx + Dy.t() * Dy) / (N - 2);
  try {
	arma::mat INV = arma::inv(COV);

	arma::mat ret = (Xxmean - Yymean) * INV * (Xxmean - Yymean).t() * nA * nU / N;
	return ret(0, 0);
  } catch (const std::runtime_error &e) {
	arma::mat INV = arma::pinv(COV);

	arma::mat ret = (Xxmean - Yymean) * INV * (Xxmean - Yymean).t() * nA * nU / N;
	return ret(0, 0);
  }
}

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

  arma::vec res = Yvec(cov.get_indices()) - cov.get_probability()(cov.get_indices());
  arma::mat ret = res.t() * K_[k] * res / 2;
  return ret(0, 0);
}

double Methods::SKATO(arma::mat &Xmat,
					  Covariates &cov,
					  arma::vec &weights,
					  const std::string &k,
					  bool shuffle,
					  int a,
					  int b) {
  // Randomize indices
  if (shuffle) {
	// arma::uvec ma_carriers = arma::find( arma::sum(Xmat) > 0 );
	cov.shuffle();
  }

  arma::uword pm = Xmat.n_cols;
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
  arma::mat Z1 = ((Xmat.each_col() % arma::sqrt(pi_1))
	  - (X1.each_col() % arma::sqrt(pi_1)) * arma::inv_sympd(X1.t() * (X1.each_col() % pi_1))
		  * (X1.t() * (Xmat.each_col() % pi_1)))
	  / arma::datum::sqrt2;

  SKAT_Optimal_GetQ skat_optimal_getQ(Xmat, res, rall);

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

double Methods::VAAST(arma::mat &Xmat,
					  arma::colvec &Yvec,
					  arma::colvec &log_casm,
					  bool score_only_minor,
					  bool score_only_alternative,
					  double site_penalty) {
  double n_case = arma::sum(Yvec);
  double n_control = Yvec.n_rows - n_case;

  arma::vec case_allele1(Xmat.n_cols);
  arma::vec control_allele1(Xmat.n_cols);

  for (arma::uword i = 0; i < Xmat.n_cols; i++) {
	case_allele1(i) = arma::dot(Yvec, Xmat.col(i));
	control_allele1(i) = arma::dot(1 - Yvec, Xmat.col(i));
  }
  arma::vec case_allele0 = 2 * n_case - case_allele1;
  arma::vec control_allele0 = 2 * n_control - control_allele1;
#if 0
  std::cerr << "case_allele0:\n" << case_allele0.t();
  std::cerr << "control_allele0:\n" << control_allele0.t();
  std::cerr << "case_allele1:\n" << case_allele1.t();
  std::cerr << "control_allele1:\n" << control_allele1.t();
#endif

  // Get ln likelihood  of each variant
  arma::vec log_lh = LRT(case_allele1, control_allele1, case_allele0, control_allele0);

  if (log_casm.n_rows == 0) {
	log_casm.zeros(Xmat.n_cols);
  }

  arma::vec vaast_site_scores = 2.0 * (log_lh + log_casm) - site_penalty;

  // mask variants with score < 1
  arma::uvec tmpmask = arma::find(vaast_site_scores <= 0);
  arma::uvec mask;
  std::vector<arma::uword> scenario;
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
	if (score_only_minor) {
	  // If both scenarios are false, append to mask
	  bool scenario1 = !(control_allele1(i) <= control_allele0(i)
		  & (case_allele1(i) * control_allele0(i) >= case_allele0(i) * control_allele1(i)));
	  bool scenario2 = !(control_allele1(i) >= control_allele0(i)
		  & (case_allele1(i) * control_allele0(i) <= case_allele0(i) * control_allele1(i)));

	  if (score_only_alternative) {
		if (scenario1)
		  scenario.push_back(i);
	  } else {
		if (scenario1 & scenario2)
		  scenario.push_back(i);
	  }
	}
  }
  mask = arma::conv_to<arma::uvec>::from(scenario);

  vaast_site_scores(mask).zeros();
#if 0
  std::cerr << "VAAST site scores:\n" << vaast_site_scores.t();
#endif
  return arma::sum(vaast_site_scores);
}

double Methods::VT(arma::mat &Xmat, arma::colvec &Yvec) {
  // All variants should be the minor allele
  arma::vec maf = arma::mean(Xmat, 0).t() / 2.;
  arma::vec hmaf = arma::unique(maf);

  arma::vec zscores(hmaf.n_rows - 1, arma::fill::zeros);
  arma::vec Ynew = Yvec - arma::mean(Yvec);

  for (int i = 0; i < hmaf.n_rows - 1; i++) {
	arma::mat Xmat_subset = Xmat.cols(arma::find(maf < hmaf[i + 1]));
	double znum = arma::sum(arma::sum(Xmat_subset.each_col() % Ynew, 1));
	double zden = std::sqrt(arma::sum(arma::sum(arma::pow(Xmat.cols(arma::find(maf < hmaf[i + 1])), 2))));
	zscores(i) = znum / zden;
  }

  return arma::max(zscores);
}

double Methods::WSS(arma::mat &Xmat, arma::colvec &Yvec) {
  double nA = arma::sum(Yvec); // Case count
  double nU = Yvec.n_rows - nA; // Control count
  double n = Yvec.n_rows;

  arma::vec mU = arma::sum(Xmat.rows(arma::find(Yvec == 0)), 0).t();
  arma::vec q = (mU + 1.) / (2. * nU + 2.);

  arma::mat w = arma::diagmat(1. / arma::sqrt(n * (q % (1. - q))));

  arma::mat gamma_mat = Xmat * w;
  gamma_mat.replace(arma::datum::nan, 0);

  arma::vec gamma = arma::sum(gamma_mat, 1);

  arma::vec ranks = rank(gamma, "ascend");

  return arma::sum(ranks(arma::find(Yvec > 0)));
}

double Methods::call(arma::mat &Xmat, arma::vec &Yvec) {
  if (method_ == "WSS") {
	return WSS(Xmat, Yvec);
  } else if (method_ == "CALPHA") {
	return CALPHA(Xmat, Yvec);
  } else if (method_ == "VT") {
	return VT(Xmat, Yvec);
  } else if (method_ == "CMC") {
	return CMC(Xmat, Yvec);
  } else if (method_ == "VAAST") {
	arma::vec log_casm;
	return VAAST(Xmat, Yvec, log_casm, false, false, 0);
  }
  return 0;
}

double Methods::call(arma::mat &Xmat, arma::vec &Yvec, arma::vec &weights) {
  if (method_ == "VAAST") {
	return VAAST(Xmat, Yvec, weights, false, false, 0);
  }
  return 0;
}

double Methods::call(arma::mat &Xmat, arma::vec &Yvec, arma::vec &weights, bool score_only_minor, double site_penalty) {
  if (method_ == "VAAST") {
	return VAAST(Xmat, Yvec, weights, score_only_minor, false, site_penalty);
  }
  throw std::logic_error("Incorrect arguments to Methods::call for method set.");
}

double Methods::call(const std::string &k, Gene &gene, Covariates &cov) {
  if (method_ == "WSS") {
	return WSS(gene.get_matrix(k), cov.get_phenotype_vector());
  } else if (method_ == "CALPHA") {
	return CALPHA(gene.get_matrix(k), cov.get_phenotype_vector());
  } else if (method_ == "VT") {
	return VT(gene.get_matrix(k), cov.get_phenotype_vector());
  } else if (method_ == "CMC") {
	return CMC(gene.get_matrix(k), cov.get_phenotype_vector(), 0);
  } else if (method_ == "VAAST") {
#if 0
	for(const auto &v : gene.get_positions(k)) {
	  std::cerr << v << "\t";
	}
	std::cerr << "\n";
#endif
	return VAAST(gene.get_matrix(k), cov.get_phenotype_vector(), gene.get_weights(k), true, true, 2);
  }
  throw (std::runtime_error("Wrong method call. 1"));
}

double Methods::call(const std::string &k, Gene &gene, Covariates &cov, arma::colvec &weights) {
  if (method_ == "WSS") {
	return WSS(gene.get_matrix(k), cov.get_phenotype_vector());
  } else if (method_ == "CALPHA") {
	return CALPHA(gene.get_matrix(k), cov.get_phenotype_vector());
  } else if (method_ == "VT") {
	return VT(gene.get_matrix(k), cov.get_phenotype_vector());
  } else if (method_ == "CMC") {
	return CMC(gene.get_matrix(k), cov.get_phenotype_vector(), 0);
  } else if (method_ == "VAAST") {
	arma::vec log_casm;
	return VAAST(gene.get_matrix(k), cov.get_phenotype_vector(), log_casm, true, true, 2);
  }
  throw (std::runtime_error("Wrong method call. 2"));
}

double Methods::call(const std::string &k, Gene &gene, Covariates &cov, bool shuffle) {
  if (method_ == "SKAT") {
	return SKAT(gene.get_matrix(k), cov, gene.get_weights(k), k, shuffle);
  } else if (method_ == "SKATO") {
	return SKATO(gene.get_matrix(k), cov, gene.get_weights(k), k, shuffle);
  }
  throw (std::runtime_error("Wrong method call. 3"));
}

std::string Methods::str() {
  return method_;
}

/*
 * VAAST Support Member Functions
 */
arma::vec Methods::LRT(arma::vec &case_allele1,
					   arma::vec &control_allele1,
					   arma::vec &case_allele0,
					   arma::vec &control_allele0) {
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

arma::vec Methods::log_likelihood(arma::vec &freq, arma::vec &allele0, arma::vec &allele1) {
  // Prevent numerical issues
  arma::vec clamped = arma::clamp(freq, 1e-9, 1.0 - 1e-9);

  return allele1 % arma::log(clamped) + allele0 % arma::log(1.0 - clamped);
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

double Methods::SKAT_Optimal_Pvalue_Davies(arma::vec &pmin_q, SKATParam &param_m, arma::vec &rall, double pmin) {
  SKAT_Integrate_Davies skat_integrate_davies(pmin_q, param_m, rall);

  // Start at DBL_EPSILON to avoid overflow
  double I;
  double error;
  try {
	auto f1 = [&](double x) { return skat_integrate_davies(x); };
	// boost::math::quadrature::tanh_sinh<double> integrator;
	// I = integrator.integrate(f1, DBL_EPSILON, 40.);
	// I = boost::math::quadrature::trapezoidal(f1, DBL_EPSILON, 40., 1e-9);
	I = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(f1, DBL_EPSILON, 40., 5, 1e-14, &error);
	// I = boost::math::quadrature::gauss<double, 20>::integrate(f1, DBL_EPSILON, 40.);
  } catch (boost::exception &e) {
	SKAT_Integrate_Liu skat_integrate_liu(pmin_q, param_m, rall);
	auto f1 = [&](double x) { return skat_integrate_liu(x); };
	// boost::math::quadrature::tanh_sinh<double> integrator;
	// I = integrator.integrate(f1, DBL_EPSILON, 40.);
	// I = boost::math::quadrature::trapezoidal(f1, DBL_EPSILON, 40., 1e-9);
	I = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(f1, DBL_EPSILON, 40., 5, 1e-14, &error);
	// I = boost::math::quadrature::gauss<double, 20>::integrate(f1, DBL_EPSILON, 40.);
  }

#ifndef NDEBUG
  std::cerr << "integral: " << I << "\n";
  std::cerr << "error: " << error << "\n";
#endif

  return 1 - I;
}

SKAT_Optimal_GetQ::SKAT_Optimal_GetQ(arma::mat &Z1, arma::vec &res, arma::vec &rall) {
  arma::uword nr = rall.n_rows;
  arma::uword pm = Z1.n_cols;

  Q_r.zeros(nr);

  arma::rowvec temp = res.t() * Z1;
  for (arma::uword i = 0; i < nr; i++) {
	double r_corr = rall(i);
	arma::mat Q1 = (1 - r_corr) * arma::sum(arma::pow(temp, 2), 1);
	arma::mat Q2 = r_corr * std::pow(pm, 2) * arma::pow(arma::mean(temp, 1), 2);
	Q_r(i) = Q1(0, 0) + Q2(0, 0);
  }
  Q_r /= 2;
}

SKATParam::SKATParam(arma::mat &Z1) {
  // Dimensions
  arma::uword n = Z1.n_rows;
  arma::uword pm = Z1.n_cols;
  arma::uword nr = 8;

  arma::vec z_mean = arma::mean(Z1, 1);
  arma::mat Z_mean = arma::mat(z_mean.n_rows, pm);
  Z_mean.each_col() = z_mean;

  arma::rowvec cof1 = (z_mean.t() * Z1) / arma::sum(arma::sum(arma::pow(z_mean, 2)));

  arma::mat Z_item1 = Z_mean * arma::diagmat(cof1);
  arma::mat Z_item2 = Z1 - Z_item1;

  arma::mat W3_2_t = Z_item2.t() * Z_item2;

  arma::vec param_lambda;
  arma::mat param_eigvec;

  arma::eig_sym(param_lambda, param_eigvec, W3_2_t);

  arma::uvec IDX1 = arma::find(param_lambda >= 0); // Positive eigenvalues
  arma::uvec IDX2 = arma::find(param_lambda > arma::mean(param_lambda(IDX1)) / 100000);

  double W3_3_item = arma::sum(arma::sum((Z_item1.t() * Z_item1) % (Z_item2.t() * Z_item2))) * 4;

  // Liu params section
  MuQ = arma::sum(param_lambda(IDX2));
  VarQ = arma::sum(arma::pow(param_lambda(IDX2), 2)) * 2 + W3_3_item;
  KerQ = arma::sum(arma::pow(param_lambda(IDX2), 4)) / std::pow(arma::sum(arma::pow(param_lambda(IDX2), 2)), 2) * 12;
  Df = 12. / KerQ;

  lambda = param_lambda;

  tau.zeros(nr);
  for (int i = 0; i < nr; i++) {
	double r_corr = rho_[i];

	double term1 = std::pow(pm, 2) * r_corr + arma::sum(arma::sum(arma::pow(cof1, 2))) * (1 - r_corr);
	tau(i) = term1 * arma::sum(arma::sum(arma::pow(z_mean, 2)));
  }
#ifndef NDEBUG
  std::cerr << "MuQ: " << MuQ << "\n";
  std::cerr << "VarQ: " << VarQ << "\n";
  std::cerr << "KerQ: " << KerQ << "\n";
  std::cerr << "Df: " << Df << "\n";
  std::cerr << "tau: " << tau.t();
#endif
}

LiuParam::LiuParam(arma::vec &c1) {
  muQ = c1(0);
  sigmaQ = std::sqrt(2 * c1(1));
  varQ = 2 * c1(1);
  s1 = c1(2) / std::pow(c1(1), 3. / 2.);
  s2 = c1(3) / std::pow(c1(1), 2);

  beta1 = std::sqrt(8) * s1;
  beta2 = 12 * s2;

  type1 = false;

  if (std::pow(s1, 2) > s2) {
	a = 1 / (s1 - std::sqrt(std::pow(s1, 2) - s2));
	d = s1 * std::pow(a, 3) - std::pow(a, 2);
	l = std::pow(a, 2) - 2 * d;
  } else {
	type1 = true;
	l = 1 / s2;
	a = std::sqrt(l);
	d = 0;
  }

  muX = l + d;
  sigmaX = std::sqrt(2) * a;
}

LiuParam::LiuParam(arma::vec &lambda, bool mod_lambda) {
  arma::vec c1{arma::sum(lambda),
			   arma::sum(arma::sum(arma::pow(lambda, 2))),
			   arma::sum(arma::sum(arma::pow(lambda, 3))),
			   arma::sum(arma::sum(arma::pow(lambda, 4)))};

  muQ = c1(0);
  sigmaQ = std::sqrt(2 * c1(1));
  varQ = 2 * c1(1);
  s1 = c1(2) / std::pow(c1(1), 3. / 2.);
  s2 = c1(3) / std::pow(c1(1), 2);

  beta1 = std::sqrt(8) * s1;
  beta2 = 12 * s2;

  type1 = false;

  if (std::pow(s1, 2) > s2) {
	a = 1 / (s1 - std::sqrt(std::pow(s1, 2) - s2));
	d = s1 * std::pow(a, 3) - std::pow(a, 2);
	l = std::pow(a, 2) - 2 * d;
  } else {
	type1 = true;
	l = 1 / s2;
	a = std::sqrt(l);
	d = 0;
  }

  muX = l + d;
  sigmaX = std::sqrt(2) * a;
}

PvalueLambda::PvalueLambda(arma::vec &lambda, double Q) {
  p_val_liu = Get_Liu_Pval_MOD_Lambda(Q, lambda);

  SKAT_Davies skd(Q, lambda);

  p_val = skd.res; // qfc in SKAT returns 1 - res;

  is_converged = true;

  // Check convergence
  if (lambda.size() == 1) {
	p_val = p_val_liu;
  } else if (skd.ifault != 0) {
	is_converged = false;
  }

  if (p_val > 1 || p_val <= 0) {
	is_converged = false;
	p_val = p_val_liu;
  }

  if (p_val == 0) {
	LiuParam param(lambda, true);
	p_val_zero_msg = Get_Liu_Pval_MOD_Lambda_Zero(Q, param);
	p_val_log = Get_Liu_Pval_MOD_Lambda(Q, lambda, true);

  }
}

double PvalueLambda::Get_Liu_Pval_MOD_Lambda(double Q, arma::vec &lambda, bool log_p) {
  LiuParam param = LiuParam(lambda, true);

  double Q_Norm = (Q - param.muQ) / param.sigmaQ;
  double Q_Norm1 = Q_Norm * param.sigmaX + param.muX;

  if (Q_Norm1 < 0) {
	if (log_p) {
	  return -arma::datum::inf;
	} else {
	  return 0;
	}
  }

  boost::math::non_central_chi_squared chisq(param.l, param.d);

  if (log_p) {
	return std::log(boost::math::cdf(chisq, Q_Norm1));
  } else {
	return boost::math::cdf(chisq, Q_Norm1);
  }
}

std::string PvalueLambda::Get_Liu_Pval_MOD_Lambda_Zero(double Q, LiuParam &param) {
  double Q_Norm = (Q - param.muQ) / param.sigmaQ;
  double Q_Norm1 = Q_Norm * param.sigmaX + param.muX;

  double temp[10] = {0.05, 1e-10, 1e-20, 1e-30, 1e-40, 1e-50, 1e-60, 1e-70, 1e-80, 1e-90};

  boost::math::non_central_chi_squared chisq(param.l, param.d);

  int max = 0;
  for (int i = 0; i < 10; i++) {
	double val = boost::math::quantile(chisq, temp[i]);
	if (val < Q_Norm1) {
	  max = i;
	}
  }
  std::stringstream ss;
  ss << "Pvalue < " << temp[max] << "\n";

  return ss.str();
}

SKAT_Davies::SKAT_Davies(double Q, arma::vec &lambda) {
  //arma::uvec rev = arma::sort_index(lambda, "descend");
  // Davies setup
  //std::vector<double> lb1 = arma::conv_to<std::vector<double>>::from(lambda(rev)); //
  std::vector<double> lb1 = arma::conv_to<std::vector<double>>::from(lambda); //
  std::vector<double> nc1(lambda.size(), 0); // nc1
  std::vector<int> df(lambda.size(), 1); // n1
  int r1 = lambda.n_rows;
  double sigma = 0; // sigma
  int lim1 = 10000;
  double acc = 10e-6;
  double trace[7]{0};
  ifault = 0;
  res = 0;


  // Davies p-value
  // c1 is Q
  // void qfc(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res){
  // qfc(&lb1[0], &nc1[0], &df[0], &r1, &sigma, &Q, &lim1, &acc, &trace[0], &ifault, &res);

  QFC qfc2(lb1, nc1, df, sigma, Q, lim1, acc);

  // res = 1 - res;
  res = 1 - qfc2.get_res();

#ifndef NDEBUG
  std::cerr << "Davies:\n";
  std::cerr << "lambda: " << lambda.t();
  std::cerr << "lb1: ";
  for(int i = 0; i < r1; i++) {
	std::cerr << lb1[i] << " ";
  }
  std::cerr << "\nnc1: ";
  for(int i = 0; i < r1; i++) {
	std::cerr << nc1[i] << " ";
  }
  std::cerr << "\ndf: ";
  for(int i = 0; i < r1; i++) {
	std::cerr << df[i] << " ";
  }
  std::cerr << "\n";
  std::cerr << "r: " << r1 << "\n";
  std::cerr << "sigma: " << sigma << "\n";
  std::cerr << "trace: ";
  for(int i = 0; i < 7; i++) {
	std::cerr << trace[i] << " ";
  }
  std::cerr << "\nifault: " << ifault << "\n";
  std::cerr << "res: " << res << "\n";
  std::cerr << "Boost approach:" << "\n";
#endif
}

SKAT_EachQ::SKAT_EachQ(arma::vec &Qall, std::vector<arma::vec> &lambda_all) {
  arma::uword nr = Qall.size();
  arma::vec c1{0, 0, 0, 0};

  std::vector<LiuParam> param_mat;

  // set vectors
  pval.zeros(nr);
  pmin_q.zeros(nr);
  pmin = 0;

  for (arma::uword i = 0; i < nr; i++) {
	double Q = Qall(i);
	arma::vec lambda_temp = lambda_all[i];

	c1(0) = arma::sum(lambda_temp);
	c1(1) = arma::sum(arma::sum(arma::pow(lambda_temp, 2)));
	c1(2) = arma::sum(arma::sum(arma::pow(lambda_temp, 3)));
	c1(3) = arma::sum(arma::sum(arma::pow(lambda_temp, 4)));

	LiuParam param_temp(c1);

	// Substitute sigmaQ for sqrt(varQ)
	double Qnorm = (Q - param_temp.muQ) / param_temp.sigmaQ * std::sqrt(2 * param_temp.l) + param_temp.l;
	PvalueLambda pval_lambda(lambda_temp, Q);

#ifndef NDEBUG
	std::cerr << "pval: " << pval_lambda.p_val << "\n";
	std::cerr << "pval_liu: " << pval_lambda.p_val_liu << "\n";
#endif

	pval(i) = pval_lambda.p_val;
	param_mat.emplace_back(param_temp);
  }

  pmin = arma::min(pval);

  for (arma::uword i = 0; i < nr; i++) {
	double muQ = param_mat[i].muQ;
	double varQ = std::pow(param_mat[i].sigmaQ, 2);
	double df = param_mat[i].l;

	boost::math::chi_squared chisq(df);

	double q_org;
	try {
	  q_org = boost::math::quantile(chisq, 1 - pmin);
	} catch (boost::exception &e) {
	  q_org = boost::math::quantile(chisq, 1 - DBL_EPSILON); // Handle overflow
	}
	double qq = (q_org - df) / std::sqrt(2 * df) * std::sqrt(varQ) + muQ;
	pmin_q(i) = qq;
  }
}

SKAT_Integrate_Davies::SKAT_Integrate_Davies(arma::vec &pmin_q, SKATParam &param_m, arma::vec &rall)
	: pmin_q(pmin_q),
	  param_m(param_m),
	  rall(rall) {}

double SKAT_Integrate_Davies::operator()(double x) {
  arma::vec temp1 = param_m.tau * x;

  arma::vec temp = (pmin_q - temp1) / (1. - rall);
  double temp_min = arma::min(temp);

  double re = 0;
  double min1 = temp_min;
  double temp_val = 0;
  if (min1 > arma::sum(param_m.lambda) * 10000) {
	temp_val = 0;
  } else {
	double min1_temp = min1 - param_m.MuQ;
	double sd1 = std::sqrt(param_m.VarQ - param_m.VarRemain) / std::sqrt(param_m.VarQ);
	double min1_st = min1_temp * sd1 + param_m.MuQ;

	SKAT_Davies skd(min1_st, param_m.lambda);

	temp_val = skd.res; // temp<-dav.re$Qq
  }
  if (temp_val > 1) {
	temp_val = 1;
  }
  boost::math::chi_squared chisq(1);
  return (1 - temp_val) * boost::math::pdf(chisq, x);
}

SKAT_Integrate_Liu::SKAT_Integrate_Liu(arma::vec &pmin_q, SKATParam &param_m, arma::vec &rall)
	: pmin_q(pmin_q),
	  param_m(param_m),
	  rall(rall) {}

double SKAT_Integrate_Liu::operator()(double x) {
  arma::mat temp1 = param_m.tau * x;

  arma::vec temp = (pmin_q - temp1) / (1 - rall);
  double temp_min = arma::min(temp);

  double temp_q = (temp_min - param_m.MuQ) / std::sqrt(param_m.VarQ) * std::sqrt(2 * param_m.Df) + param_m.Df;

  boost::math::chi_squared chisq(param_m.Df);
  boost::math::chi_squared chisq2(1);
  return boost::math::cdf(chisq, temp_q) * boost::math::pdf(chisq2, x);
}
