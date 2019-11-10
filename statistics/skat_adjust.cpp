//
// Created by Bohlender,Ryan James on 8/23/18.
//

#include "skat_adjust.hpp"

#include <memory>

// Boost Math
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/chi_squared.hpp>

// Boost Integration
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>

constexpr std::array<double, 8> SKAT_Optimal_Logistic_VarMatching::rcorr;

SKAT_Residuals_Logistic::SKAT_Residuals_Logistic(arma::mat &cov, arma::vec &mu, arma::vec &Y, int nresampling) {
  nResampling = static_cast<arma::uword>(nresampling);
  arma::mat covt = cov.t();
  calcX1(covt);

  res = Y - mu;
  if(nResampling > 0) {
    res_out.reshape(res.n_rows, nResampling);
    for(arma::uword i = 0; i < nResampling; i++) {
      res_out.col(i) = arma::shuffle(res);
    }
  }

  pi_1 = arma::diagmat(mu) * (1 - mu);
  this->mu = mu;
}

void SKAT_Residuals_Logistic::calcX1(arma::mat &cov) {
  arma::mat Q, R;

  arma::qr(Q, R, cov);

  if (arma::rank(Q) <= cov.n_cols) {
	arma::mat U, V;
	arma::vec s;

	svd_econ(U, s, V, cov, "left");

	X1 = U;
  } else {
	X1 = cov;
  }
}

SKAT_Residuals_Logistic::SKAT_Residuals_Logistic(Covariates &cov, int nresampling) {
  nResampling = static_cast<arma::uword>(nresampling);
  arma::mat covt = cov.get_covariate_matrix().t();
  calcX1(covt);

  res = cov.get_residuals();

  if(nResampling > 0) {
	res_out.reshape(res.n_rows, nResampling);
	for(arma::uword i = 0; i < nResampling; i++) {
	  res_out.col(i) = arma::shuffle(res);
	}
  }

  mu = cov.get_fitted();
  pi_1 = arma::diagmat(mu) * (1 - mu);
}

SKAT_Adjust::SKAT_Adjust(Gene &gene,
						 Covariates &cov,
						 const std::string &k,
						 const std::string &kernel,
						 int a,
						 int b,
						 std::shared_ptr<SKAT_Residuals_Logistic> &re2_,
						 std::map<std::string, std::shared_ptr<SKAT_Optimal_GetQ>> &Q_sim_all) {
  arma::sp_mat Z(gene.get_matrix(k));

  if (gene.get_weights(k).n_rows == 0) {
	weights.reshape(Z.n_cols, 1);
	arma::vec maf(arma::mean(Z, 0).t() / 2.);

	for (arma::uword i = 0; i < Z.n_cols; i++) {
	  weights(i) = std::pow(maf(i), a - 1) * std::pow(1 - maf(i), b - 1) / boost::math::beta(a, b);
	}
  }

  arma::vec phen_vec = cov.get_phenotype_vector();
  arma::vec prob_vec = cov.get_fitted();

  // Residuals
  re1 = SKAT_Residuals_Logistic(cov.get_covariate_matrix(), prob_vec, phen_vec, 0);
  if(re2_ == nullptr) {
	re2_ = std::make_shared<SKAT_Residuals_Logistic>(cov.get_covariate_matrix(), prob_vec, phen_vec, 10000);
  }

  SKAT_Optimal_Logistic_VarMatching re(re1,
									   *re2_,
									   Q_sim_all[k],
									   Z,
									   kernel,
									   weights);

  p_value = re.pval;
}

/**
 * @brief
 * @param re1 The residuals.
 * @param re2 The residuals and permuted residuals for moment adjustment.
 * @param Q_sim_all Q from resamples for moment adjustment.
 * @param Z Genotype matrix
 * @param kernel Which kernel to use.
 * @param weights Optional weight vector.
 */
SKAT_Optimal_Logistic_VarMatching::SKAT_Optimal_Logistic_VarMatching(SKAT_Residuals_Logistic &re1,
																	 SKAT_Residuals_Logistic &re2,
																	 std::shared_ptr<SKAT_Optimal_GetQ> &Q_sim_all,
																	 arma::sp_mat &Z,
																	 const std::string &kernel,
																	 arma::vec &weights) {
  arma::uword n = Z.n_rows;
  arma::uword pm = Z.n_cols;
  arma::uword nr = rcorr.size();

  arma::vec rall(nr);
  for (arma::uword i = 0; i < nr; i++) {
	rall(i) = rcorr[i];
  }

  arma::vec &pi_1 = re1.pi_1; // References to shorten calls;
  arma::mat X1 = re1.X1;     // References to shorten calls;

  // TODO Add checks for other kernels and throw error
  if (kernel == "wLinear") {
	Z = arma::diagmat(weights) * Z;
  }

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

  SKAT_Optimal_GetQ outQ(Z, re1.res, rall, re1.res_out, re1.nResampling); // Normal Q calculation
  if(Q_sim_all == nullptr) {
	Q_sim_all = std::make_shared<SKAT_Optimal_GetQ>(Z, re1.res, rall, re2.res_out, re2.nResampling); // For moment adjustment
  }

  // P-Values
  arma::vec p_all = re1.mu;

  SKAT_Optimal_Get_Pvalue_VarMatching out(outQ.Q_r,
										  Z1,
										  rall,
										  p_all,
										  (*Q_sim_all).Q_sim,
										  re1.res,
										  re2.res_out);

  pval_each = out.pval_each;
  qval_each = outQ.Q_r;
  minp = arma::min(pval_each);

  arma::uvec id_temp = arma::find(pval_each == minp);
  arma::uvec id_temp1 = arma::find(rall >= 0.999);

  rho_est = rcorr[id_temp(0)];

  pval = out.pval;
}

SKAT_Optimal_Get_Pvalue_VarMatching::SKAT_Optimal_Get_Pvalue_VarMatching(arma::vec &Q_all,
																		 arma::mat &Z1,
																		 arma::vec &rall,
																		 arma::vec &p_all,
																		 arma::mat &Q_sim_all,
																		 arma::vec &res,
																		 arma::mat &res_moments) {
  arma::uword nr = rall.n_rows;
  arma::uword nq = 1; // Number of resamples -- we resample elsewhere
  arma::uword pm = Z1.n_cols;

  // lambda.all
  // Z2.all
  std::vector<arma::mat> Z2_all;
  for (arma::uword i = 0; i < nr; i++) {
	double rcorr = rall(i);

	arma::mat Rm = arma::diagmat(arma::vec(pm).fill(1 - rcorr)) + arma::mat(pm, pm).fill(rcorr);
	Z2_all.emplace_back(Z1 * arma::chol(Rm).t());
  }

  // Mixture parameters
  param_m = SKAT_Optimal_Param_VarMatching(Z1, rall, p_all, res, res_moments);
  each_info = SKAT_Optimal_Each_Q_VarMatching(param_m,
											  Q_all,
											  rall,
											  Z2_all,
											  p_all,
											  Q_sim_all);

  arma::vec pmin_q = each_info.pmin_q;
  double pmin = each_info.pmin;

  double muQ = param_m.param.muQ;
  double varQ = param_m.param.varQ + param_m.VarRemain;
  double df = param_m.param.df;
  arma::vec tau = param_m.tau;

  pval_each = each_info.pval;

  pval = SKAT_Optimal_Pvalue_VarMatching(pmin_q, muQ, varQ, df, tau, rall, pmin);

  double multi = 3;

  arma::uvec IDX = arma::find(pval_each > 0);

  double pval1 = arma::min(pval_each) * multi;
  if (pval < 0 || IDX.n_rows < rall.n_rows) {
	pval = pval1;
  }

  if (pval == 0) {
	if (IDX.n_rows > 0) {
	  pval = arma::min(pval_each(IDX));
	}
  }
}

double SKAT_Optimal_Get_Pvalue_VarMatching::SKAT_Optimal_Pvalue_VarMatching(arma::vec pmin_q,
																			double muQ,
																			double varQ,
																			double df,
																			arma::vec tau,
																			arma::vec &rall,
																			double pmin) {
  // SKAT_Optimal_Integrate_Func_VarMatching
  auto f1 = [&](double x) {
	arma::uword nr = rall.size();

	arma::vec temp1 = tau * x;
	arma::vec temp = (pmin_q - temp1) / (1. - rall);

	double temp_min = arma::min(temp);
	double temp_q, val;

	if (varQ > 0) {
	  temp_q = (temp_min - muQ) / std::sqrt(varQ) * std::sqrt(2. * df) + df;

	  boost::math::chi_squared chisq(df);
	  try {
		val = boost::math::cdf(chisq, temp_q);
	  } catch (boost::exception &e) {
		// input < 0
		val = 0;
	  }
	} else {
	  val = 0;
	}
	boost::math::chi_squared chisq(1);
	double re = val * boost::math::pdf(chisq, x);
	return re;
  };

  double error;
  double I = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(f1, DBL_EPSILON, 40., 5, 1e-14, &error);
  // boost::math::quadrature::tanh_sinh<double> integrator;
  // double I = integrator.integrate(f1, DBL_EPSILON, 40.);

  double pvalue = 1. - I;

  if (!std::isnan(pmin)) {
	if (pmin * rall.size() < pvalue) {
	  pvalue = pmin * rall.size();
	}
  }
  return pvalue;

}

SKAT_Optimal_Param_VarMatching::SKAT_Optimal_Param_VarMatching(arma::mat &Z1,
															   arma::vec &rall,
															   arma::vec &p_all,
															   arma::vec &res,
															   arma::mat &res_moments) {
  arma::uword n = Z1.n_rows;
  arma::uword pm = Z1.n_cols;
  arma::uword nr = rall.n_rows;

  arma::vec z_mean = arma::mean(Z1, 1);
  arma::mat Z_mean = arma::mat(z_mean.n_rows, pm);
  Z_mean.each_col() = z_mean;

  arma::rowvec cof1 = (z_mean.t() * Z1) / arma::sum(arma::sum(arma::pow(z_mean, 2)));

  arma::mat Z_item1 = Z_mean * arma::diagmat(cof1);
  arma::mat Z_item2 = Z1 - Z_item1;

  arma::mat Q_temp_res1 = res_moments.t() * Z_item2;
  Q_sim = arma::sum(arma::pow(Q_temp_res1, 2), 1) / 2.;

  param = SKAT_Logistic_VarMatching_Param(Z_item2, p_all, Q_sim);

  // W3.3.item
  VarRemain = arma::sum(arma::sum((Z_item1.t() * Z_item1) % (Z_item2.t() * Z_item2))) * 4;

  tau.reshape(nr, 1);
  for (arma::uword i = 0; i < nr; i++) {
	double r_corr = rall(i);
	arma::rowvec term1 = pm * r_corr + arma::pow(cof1, 2) * (1 - r_corr);
	tau(i) = arma::sum(term1) * arma::sum(arma::sum(arma::pow(z_mean, 2)));
  }
}

SKAT_Optimal_Param_VarMatching &SKAT_Optimal_Param_VarMatching::operator=(const SKAT_Optimal_Param_VarMatching &rhs) {
  Q_sim = rhs.Q_sim;
  VarRemain = rhs.VarRemain;
  tau = rhs.tau;
  param = rhs.param;

  return *this;
}

SKAT_Logistic_VarMatching_Param::SKAT_Logistic_VarMatching_Param(arma::mat &Z1, arma::vec &p_all, arma::vec &Q_sim)
: muQ(0),
  varQ(0),
  df(0),
  param_noadj(),
  nlambda(0) {
  //	re.param<-SKAT_Logistic_VarMatching_GetParam1(Z.item2, p_all, Q.sim, type)

  auto lambda_u = Get_Lambda_U_From_Z(Z1);

  if (!lambda_u.first.empty()) {
	calc_param(lambda_u.first, lambda_u.second, p_all, Q_sim);
  } else {
	only_sim(Z1, p_all, Q_sim);
  }
  // End re.param<-SKAT_Logistic_VarMatching_GetParam1(Z.item2, p_all, Q.sim, type)
}

/**
 * @brief Calculates svd on Z1
 * @param Z1 Matrix
 * @return Pair containing .first = lambda, .second = U
 */
std::pair<arma::vec, arma::mat> SKAT_Logistic_VarMatching_Param::Get_Lambda_U_From_Z(arma::mat &Z1) {
  arma::uword pm = Z1.n_cols;
  std::pair<arma::mat, arma::mat> ret;

  if (pm == 1) {
	ret.first = arma::sum(arma::sum(arma::pow(Z1, 2)));
	ret.second = Z1 / arma::sqrt(ret.first);
	return ret;
  }

  arma::mat U, V;
  arma::vec s; // Equal to d in R's svd
  arma::svd(U, s, V, Z1);

  arma::vec lambda_org = arma::pow(s, 2);

  arma::uvec IDX = arma::find(lambda_org > arma::mean(lambda_org) / 100000.);

  if (IDX.size() == 0) {
	return ret;
  }

  ret.first = lambda_org(IDX);
  ret.second = U.cols(IDX);

  return ret;
}

/**
 * @brief Calculate variance entry in cov matrix
 * @param m4 m4 vector
 * @param p_all mu
 * @param Ui i'th eigenvector
 * @param Uj j'th eigenvector
 * @return double -- variance entry
 */
double SKAT_Logistic_VarMatching_Param::SKAT_Get_Var_Elements(arma::vec &m4,
															  arma::vec &p_all,
															  arma::subview_col<double> Ui,
															  arma::subview_col<double> Uj) {
  arma::vec temp1 = arma::diagmat(arma::pow(Ui, 2)) * arma::pow(Uj, 2);

  double a1 = arma::sum(arma::sum(arma::diagmat(m4) * temp1));
  double a2 = arma::sum(arma::sum(arma::pow(Ui, 2))) * arma::sum(arma::sum(arma::pow(Uj, 2))) - arma::sum(temp1);
  double a3 = (std::pow(arma::sum(arma::sum(arma::diagmat(Ui) * Uj)), 2) - sum(temp1)) * 2;

  return a1 + a2 + a3;
}

double SKAT_Logistic_VarMatching_Param::SKAT_Get_DF_Sim(arma::vec &Q_sim) {
  double s2_sim = SKAT_Get_Kurtosis(Q_sim);
  double df_sim = 12 / s2_sim;

  if (s2_sim <= 0) {
	df_sim = 100000;
  } else if (df_sim < 0.01) {
	double s1_sim = SKAT_Get_Skewness(Q_sim);
	df_sim = 8. / std::pow(s1_sim, 2);
  }

  return df_sim;
}

double SKAT_Logistic_VarMatching_Param::SKAT_Get_Kurtosis(arma::vec &x) {
  if (arma::stddev(x) == 0) {
	return -100;
  }

  double m4 = arma::mean(arma::pow(x - arma::mean(x), 4));
  return m4 / std::pow(arma::stddev(x), 4) - 3;
}

double SKAT_Logistic_VarMatching_Param::SKAT_Get_Skewness(arma::vec &x) {
  double m3 = arma::mean(arma::pow(x - arma::mean(x), 3));

  return m3 / std::pow(arma::stddev(x), 3);
}

/**
 * @brief Get parameters for p-value calculation - SKAT_Get_Cov_Param
 * @param lambda eigenvalues
 * @param U eigenvectors
 * @param p_all mu
 * @param Q_sim
 */
void SKAT_Logistic_VarMatching_Param::calc_param(arma::vec &lambda, arma::mat &U, arma::vec &p_all, arma::vec &Q_sim) {
  arma::uword pm = lambda.size();
  arma::vec p_all_prod = arma::diagmat(p_all) * (1 - p_all);
  arma::vec m4 = arma::diagmat(p_all_prod) * (3 * arma::pow(p_all, 2) - 3 * p_all + 1) / arma::pow(p_all_prod, 2);

  // SKAT_Get_Cov_Param
  zeta = arma::vec(pm, arma::fill::zeros);
  var_i = arma::vec(pm, arma::fill::zeros);

  for (arma::uword i = 0; i < pm; i++) {
	double temp_M1 =
		std::pow(arma::sum(arma::sum(arma::pow(U.col(i), 2))), 2) - arma::sum(arma::sum(arma::pow(U.col(i), 4)));
	zeta(i) = arma::sum(arma::sum(arma::diagmat(m4) * arma::pow(U.col(i), 4))) + 3 * temp_M1;
	var_i(i) = zeta(i) - 1;
  }

  arma::mat cov_mat;
  if (pm == 1) {
	cov_mat.reshape(1, 1);
	cov_mat = zeta * arma::pow(lambda, 2);
  } else {
	cov_mat = arma::diagmat(arma::diagmat(zeta) * arma::pow(lambda, 2));
	for (arma::uword i = 0; i < pm - 1; i++) {
	  for (arma::uword j = i + 1; j < pm; j++) {
		cov_mat(i, j) = SKAT_Get_Var_Elements(m4, p_all, U.col(i), U.col(j));
		cov_mat(i, j) = cov_mat(i, j) * lambda(i) * lambda(j);
	  }
	}
  }

  cov_mat = cov_mat + cov_mat.t();
  cov_mat.diag() /= 2;

  muQ = arma::sum(lambda);
  varQ = arma::sum(arma::sum(cov_mat)) - std::pow(arma::sum(lambda), 2);
  lambda_new = arma::diagmat(lambda) * arma::sqrt(var_i) / arma::datum::sqrt2;
  nlambda = lambda_new.n_rows;
  // End SKAT_Get_Cov_Param

  double
	  s2 = arma::sum(arma::sum(arma::pow(lambda_new, 4))) / std::pow(arma::sum(arma::sum(arma::pow(lambda_new, 2))), 2);
  df = 1. / s2;

  if (Q_sim.n_rows > 0) {
	df = SKAT_Get_DF_Sim(Q_sim);
  }

  arma::vec c1(4);
  for (arma::uword i = 0; i < 4; i++) {
	c1(i) = arma::sum(arma::sum(arma::pow(lambda, i)));
  }

  param_noadj = LiuParam(c1);
}

void SKAT_Logistic_VarMatching_Param::only_sim(arma::mat &Z1, arma::vec &p_all, arma::vec &Q_sim) {
  arma::mat tempZ1 = Z1.t() * Z1;
  arma::vec lambda = Get_Lambda(tempZ1);

  // No Adjustment
  arma::vec c1(4, arma::fill::zeros);
  for (arma::uword i = 0; i < 4; i++) {
	c1(i) = arma::sum(arma::sum(arma::pow(lambda, i)));
  }

  muQ = arma::sum(lambda);
  varQ = arma::var(Q_sim);
  df = SKAT_Get_DF_Sim(Q_sim);

  param_noadj = LiuParam(c1);
}

SKAT_Optimal_Each_Q_VarMatching::SKAT_Optimal_Each_Q_VarMatching(SKAT_Optimal_Param_VarMatching &param_m,
																 arma::vec &Q_all,
																 arma::vec &rall,
																 std::vector<arma::mat> &Z2_all,
																 arma::vec &p_all,
																 arma::mat &Q_sim_all)
																 : pmin(0) {

  arma::uword nr = rall.size();
  arma::vec c1(4, arma::fill::zeros);

  pval.zeros(nr);
  pmin_q.zeros(nr);

  std::vector<SKAT_Logistic_VarMatching_Param> re_param;
  for (arma::uword i = 0; i < nr; i++) {
	double Q = Q_all(i);

	SKAT_PValue_Logistic_VarMatching out(Q, Z2_all[i], p_all, Q_sim_all.col(i));

	re_param.emplace_back(out.param);
	pval(i) = out.p_value;
  }

  pmin = arma::min(pval);
  // Readjust kurtosis of Q
  for (arma::uword i = 0; i < nr; i++) {
	double r_corr = rall(i);
	double muQ = re_param[i].muQ;
	double varQ = std::pow(1 - r_corr, 2) * (param_m.param.varQ + param_m.VarRemain) + std::pow(param_m.tau(i), 2) * 2.;

	double vq1 = param_m.param.varQ + param_m.VarRemain;
	double kur = SKAT_Optimal_Kurtosis_Mixture(param_m.param.df, 1, vq1, (1 - r_corr), param_m.tau[i]);

	double df = 12. / kur;

	boost::math::chi_squared chisq(df);
	double q_org = boost::math::quantile(chisq, std::max(1. - pmin, DBL_EPSILON));
	double qq = (q_org - df) / std::sqrt(2 * df) * std::sqrt(varQ) + muQ;

	pmin_q(i) = qq;
  }
}

double SKAT_Optimal_Each_Q_VarMatching::SKAT_Optimal_Kurtosis_Mixture(double df1,
																	  double df2,
																	  double v1,
																	  double a1,
																	  double a2) {
  double v2 = 2 * df2;

  double S4_1 = (12. / df1 + 3.) * std::pow(v1, 2);
  double S4_2 = (12. / df2 + 3.) * std::pow(v2, 2);

  double S4 = std::pow(a1, 4) * S4_1 + std::pow(a2, 4) * S4_2 + 6 * std::pow(a1, 2) * std::pow(a2, 2) * v1 * v2;
  double S2 = std::pow(a1, 2) * v1 + std::pow(a2, 2) * v2;

  double K = S4 / std::pow(S2, 2) - 3.;

  if (K < 0) {
	K = 0.0001;
  }

  return K;
}

SKAT_PValue_Logistic_VarMatching::SKAT_PValue_Logistic_VarMatching(double Q,
																   arma::mat &Z1,
																   arma::vec &p_all,
																   arma::subview_col<double> Q_sim) {
  arma::vec Qsim = Q_sim; // Get concrete type
  param = SKAT_Logistic_VarMatching_Param(Z1, p_all, Qsim);
  std::string pval_msg;

  if (param.varQ == 0) {
	p_value = 1;
  } else {
	boost::math::non_central_chi_squared chisq(param.df, 0);
	double Q_norm = (Q - param.muQ) / std::sqrt(param.varQ);
	double Q_norm1 = Q_norm * std::sqrt(2 * param.df) + param.df;
	try {
	  p_value = 1. - boost::math::cdf(chisq, Q_norm1);
	} catch (boost::exception &e) {
	  p_value = 1;
	}
  }

  if (param.param_noadj.sigmaQ == 0) {
	p_value_noadj = 1;
  } else {
	double Q_norm = (Q - param.param_noadj.muQ) / param.param_noadj.sigmaQ;
	double Q_norm1 = Q_norm * param.param_noadj.sigmaX + param.param_noadj.muX;

	boost::math::non_central_chi_squared chisq(param.param_noadj.l, 0);
	try {
	  p_value_noadj = 1 - boost::math::cdf(chisq, Q_norm1);
	} catch (boost::exception &e) {
	  p_value_noadj = 1;
	}
  }
}

arma::vec Get_Lambda(arma::mat &K) {
  // eigvec is ignored
  arma::vec lambda1;
  arma::mat eigvec;

  arma::eig_sym(lambda1, eigvec, K);

  arma::uvec IDX1 = arma::find(lambda1 >= 0);
  arma::uvec IDX2 = arma::find(lambda1 > arma::mean(lambda1(IDX1)) / 100000);

  return lambda1(IDX2);
}
