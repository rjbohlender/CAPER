//
// Created by Bohlender,Ryan James on 8/22/18.
//

#include "skat.hpp"
#include "../third_party/QFC/qfc2.hpp"

// Boost Math
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/chi_squared.hpp>

// Boost Integration
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>

constexpr double SKATParam::rho_[];

SKAT_Optimal_GetQ::SKAT_Optimal_GetQ(arma::sp_mat &Z,
									 arma::vec &res,
									 arma::vec &rall,
									 arma::mat &res_out,
									 int nResampling) {
  arma::uword nr = rall.n_rows;
  arma::uword pm = Z.n_cols;

  Q_r.zeros(nr);

  arma::rowvec temp = res.t() * Z;
  for (arma::uword i = 0; i < nr; i++) {
	double r_corr = rall(i);
	arma::mat Q1 = (1 - r_corr) * arma::sum(arma::pow(temp, 2), 1);
	arma::mat Q2 = r_corr * std::pow(pm, 2) * arma::pow(arma::mean(temp, 1), 2);
	Q_r(i) = Q1(0, 0) + Q2(0, 0);
  }
  Q_r = Q_r / 2.;

  // For moment adjustment
  if (nResampling > 0) {
	Q_sim = arma::mat(nResampling, nr);

	for (arma::uword i = 0; i < nResampling; i++) {
	  if (res_out.size() > 0) {
		temp = res_out.col(i).t() * Z;
	  } else {
		temp = arma::shuffle(res).t() * Z;
	  }
	  for (arma::uword j = 0; j < nr; j++) {
		double r_corr = rall(j);
		arma::mat Q1 = (1 - r_corr) * arma::sum(arma::pow(temp, 2), 1);
		arma::mat Q2 = r_corr * std::pow(pm, 2) * arma::pow(arma::mean(temp, 1), 2);
		Q_sim(i, j) = Q1(0, 0) + Q2(0, 0);
	  }
	}
	Q_sim = Q_sim / 2.;
  }
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

  // Make Z_item2 symmetric
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
  for (arma::uword i = 0; i < nr; i++) {
	double r_corr = rho_[i];

	double term1 = std::pow(pm, 2) * r_corr + arma::sum(arma::sum(arma::pow(cof1, 2))) * (1 - r_corr);
	tau(i) = term1 * arma::sum(arma::sum(arma::pow(z_mean, 2)));
  }
#if 0
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
  for (arma::uword i = 0; i < 10; i++) {
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
  int lim1 = 1000000;
  double acc = 10e-9;
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

#if 0
  std::cerr << "Davies:\n";
  std::cerr << "lambda: " << lambda.t();
  std::cerr << "lb1: ";
  for (int i = 0; i < r1; i++) {
	std::cerr << lb1[i] << " ";
  }
  std::cerr << "\nnc1: ";
  for (int i = 0; i < r1; i++) {
	std::cerr << nc1[i] << " ";
  }
  std::cerr << "\ndf: ";
  for (int i = 0; i < r1; i++) {
	std::cerr << df[i] << " ";
  }
  std::cerr << "\n";
  std::cerr << "r: " << r1 << "\n";
  std::cerr << "sigma: " << sigma << "\n";
  std::cerr << "trace: ";
  for (int i = 0; i < 7; i++) {
	std::cerr << trace[i] << " ";
  }
  std::cerr << "\nifault: " << ifault << "\n";
  std::cerr << "res: " << res << "\n";
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

#if 0
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

double SKAT_Optimal_Pvalue_Davies(arma::vec &pmin_q, SKATParam &param_m, arma::vec &rall, double pmin) {
  SKAT_Integrate_Davies skat_integrate_davies(pmin_q, param_m, rall);

  // Start at DBL_EPSILON to avoid overflow
  double I;
  double error;
  try {
	auto f1 = [&](double x) { return skat_integrate_davies(x); };
	// boost::math::quadrature::tanh_sinh<double> integrator;
	// I = integrator.integrate(f1, DBL_EPSILON, 40.);
	// I = boost::math::quadrature::trapezoidal(f1, DBL_EPSILON, 40., 1e-9);
	I = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(f1, DBL_EPSILON, 40., 15, 1e-14, &error);
	// I = boost::math::quadrature::gauss<double, 20>::integrate(f1, DBL_EPSILON, 40.);
  } catch (boost::exception &e) {
	SKAT_Integrate_Liu skat_integrate_liu(pmin_q, param_m, rall);
	auto f1 = [&](double x) { return skat_integrate_liu(x); };
	// boost::math::quadrature::tanh_sinh<double> integrator;
	// I = integrator.integrate(f1, DBL_EPSILON, 40.);
	// I = boost::math::quadrature::trapezoidal(f1, DBL_EPSILON, 40., 1e-9);
	I = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(f1, DBL_EPSILON, 40., 15, 1e-14, &error);
	// I = boost::math::quadrature::gauss<double, 20>::integrate(f1, DBL_EPSILON, 40.);
  }

#if 0
  std::cerr << "integral: " << I << "\n";
  std::cerr << "error: " << error << "\n";
#endif

  return 1 - I;
}

