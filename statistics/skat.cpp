//
// Created by Bohlender,Ryan James on 9/22/18.
//

#include "skat.hpp"

#include "../third_party/QFC/qfc2.hpp"

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>

#include <boost/math/tools/roots.hpp>
#include <utility>

double SKAT_pval(double Q, const arma::vec &lambda, bool sadd) {
  if (std::isnan(Q)) {
    return 1;
  }

  if (!sadd) {
    std::vector<double> lb1 = arma::conv_to<std::vector<double>>::from(lambda);
    std::vector<double> nc1(lb1.size(), 0);
    std::vector<int> df(lb1.size(), 1);
    double sigma = 0;
    int lim1 = 1000000;
    double acc = 1e-9;

    QFC qfc(lb1, nc1, df, sigma, Q, lim1, acc);

    double pval = 1. - qfc.get_res();
    if (pval >= 1 || pval <= 0 || qfc.get_fault() > 0) {
      pval = Saddlepoint(Q, lambda);
    }
    return pval;
  } else {
    double pval = Saddlepoint(Q, lambda);
    return pval;
  }
}

double Liu_qval_mod(double pval, const arma::vec &lambda) {
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

  double a, d, l;
  if (std::pow(s1, 2) > s2) {
    a = 1. / (s1 - std::sqrt(s1 * s1 - s2));
    d = s1 * std::pow(a, 3) - std::pow(a, 2);
    l = std::pow(a, 2) - 2 * d;
  } else { // Modified to match kurtosis only in this branch
    l = 1. / s2;
    a = std::sqrt(l);
    d = 0;
  }

  double muX = l + d;
  double sigmaX = arma::datum::sqrt2 * a;
  double df = l;

  boost::math::chi_squared chisq(df);
  double q = boost::math::quantile(boost::math::complement(chisq,
                                                           pval > 0 ? pval
                                                                    : std::sqrt(std::numeric_limits<double>::min())));
  return (q - muX) / sigmaX * sigmaQ + muQ; // Does match the Liu (2009) paper.
}


double Liu_pval(double Q, const arma::vec &lambda) {
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

  double a, d, l; // l = degrees of freedom, d = non-centrality parameter for non-central chisq
  if (std::pow(s1, 2) > s2) {
    a = 1. / (s1 - std::sqrt(std::pow(s1, 2) - s2));
    d = s1 * std::pow(a, 3) - std::pow(a, 2);
    l = std::pow(a, 2) - 2 * d;
  } else {
    a = std::sqrt(1. / s2);
    d = 0;
    l = 1. / s2;
  }
  double muX = l + d;
  double sigmaX = arma::datum::sqrt2 * a;
  double df = l;

  double Qnorm = (Q - muQ) / sigmaQ * sigmaX + muX;

  try {
    boost::math::non_central_chi_squared chisq(df, d);
    return boost::math::cdf(boost::math::complement(chisq, Qnorm));
  } catch (std::exception &e) {
    return 1;
  }
}

double Saddlepoint(double Q, const arma::vec &lambda) {
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
  auto kppprime0 = [&](double &zeta) -> double {
    return 8 * arma::accu(arma::pow(ulambda, 3) / arma::pow(1 - 2 * (zeta * ulambda), 3));
  };
#if 0
  auto hatzetafn = [&](double zeta) -> double {
    return kprime0(zeta) - Q;
  };
#else
  auto hatzetafn = [&](double zeta) -> std::tuple<double, double, double> {
    return std::make_tuple(kprime0(zeta) - Q, kpprime0(zeta), kppprime0(zeta));
  };
#endif

  arma::uword n = ulambda.size();

  double lmin, lmax;
  if (arma::any(ulambda < 0)) {
    lmin = arma::max(1 / (2 * ulambda(arma::find(ulambda < 0)))) * 0.99999;
  } else if (Q > arma::sum(ulambda)) {
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
      tmp = boost::math::tools::bisect(hatzetafn, lmin, lmax, tol, max_iter);
#elif 0
  int digits = std::numeric_limits<double>::digits - 3;
  boost::math::tools::eps_tolerance<double> tol(digits);
  boost::uintmax_t max_iter = 1000;
  double factor = (lmax - lmin) / 100;
  std::pair<double, double>
      tmp = boost::math::tools::bracket_and_solve_root(hatzetafn, lmin, factor, true, tol, max_iter);
#else
  int digits = std::numeric_limits<double>::digits - 3;
  // boost::math::tools::eps_tolerance<double> tol(digits);
  int tol = static_cast<int>(digits * 0.6);
  boost::uintmax_t max_iter = 1000;
  double guess = (lmax - lmin) / 2.;
  try {
#ifdef __clang__
    double hatzeta = boost::math::tools::halley_iterate(hatzetafn, guess, lmin, lmax, tol, max_iter);
#else
    double hatzeta = boost::math::tools::schroder_iterate(hatzetafn, guess, lmin, lmax, tol, max_iter);
#endif
    // double hatzeta = tmp.first + (tmp.second - tmp.first) / 2.;

    double w = sgn(hatzeta) * std::sqrt(2 * (hatzeta * Q - k0(hatzeta)));
    double v = hatzeta * std::sqrt(kpprime0(hatzeta));

    if (std::abs(hatzeta) < 1e-4 || std::isnan(w) || std::isnan(v)) {
      return Liu_pval(Q * d, lambda);
    } else {
      boost::math::normal norm(0, 1);
      return boost::math::cdf(boost::math::complement(norm, w + std::log(v / w) / w));
    }
  } catch (boost::math::evaluation_error &e) {
    return SKAT_pval(Q, lambda, false);
  }
#endif

}


SKATR_Null::SKATR_Null(Covariates cov)
: crand(std::random_device{}()), cov_(std::move(cov)) {
  X = cov_.get_covariate_matrix();
  Y = cov_.get_phenotype_vector();
  pi0 = cov_.get_fitted();
  Yv = pi0 % (1 - pi0);
  Yh = arma::sqrt(Yv);

  arma::mat U, V;
  arma::vec S;
  arma::svd_econ(U, S, V, arma::diagmat(Yh) * X);

  Ux = arma::diagmat(Yh) * U;
  U0 = Y - pi0;
}

auto SKATR_Null::shuffle(arma::vec &phenotypes) -> void {
  cov_.set_phenotype_vector(phenotypes);
  cov_.refit_permuted();

  pi0 = cov_.get_fitted();
  Yv = pi0 % (1 - pi0);
  Yh = arma::sqrt(Yv);
  Y = phenotypes;

  arma::mat U, V;
  arma::vec S;
  bool success = arma::svd_econ(U, S, V, arma::diagmat(Yh) * X, "left");
  if (!success) {
    std::cerr << "Yh: " << Yh.t();
    std::cerr << "Yv: " << Yv.t();
    std::cerr << "Y: " << Y.t();
    for (int i = 0; i < Y.n_elem; i++) {
      std::cerr << Y[i];
      for (int j = 0; j < X.n_cols; j++) {
        std::cerr << "\t" << X(i, j);
      }
      std::cerr << std::endl;
    }
  }

  Ux = arma::diagmat(Yh) * U;
  U0 = Y - pi0;
}

auto SKATR_Null::get_U0() noexcept -> arma::vec {
  return U0;
}

auto SKATR_Null::get_pi0() noexcept -> arma::vec {
  return pi0;
}

auto SKATR_Null::get_Yv() noexcept -> arma::vec {
  return Yv;
}

auto SKATR_Null::get_Yh() noexcept -> arma::vec {
  return Yh;
}

auto SKATR_Null::get_Ux() noexcept -> const arma::mat & {
  return Ux;
}

auto SKATR_Null::get_coef() noexcept -> arma::rowvec {
  return coef;
}

auto SKATR_Null::get_Y() noexcept -> arma::vec {
  return Y;
}

auto SKATR_Null::get_X() noexcept -> const arma::mat & {
  return X;
}

SKATR_Linear_Null::SKATR_Linear_Null(Covariates cov)
: cov_(std::move(cov)) {
  X = cov_.get_covariate_matrix();
  Y = cov_.get_phenotype_vector();

  arma::mat U, V;
  arma::vec s;
  arma::svd_econ(U, s, V, X.t());

  Ux = U;
  U0 = cov_.get_residuals();
  s2 = arma::accu(arma::pow(U0, 2)) / (Y.n_elem - X.n_rows); // Residual sum of squares over residual df
}

auto SKATR_Linear_Null::shuffle(arma::vec &phenotypes) -> void {
  cov_.set_phenotype_vector(phenotypes);
  cov_.refit_permuted();

  Y = phenotypes;

  U0 = cov_.get_residuals();
  s2 = arma::accu(arma::pow(U0, 2)) / (Y.n_elem - X.n_rows);
}

auto SKATR_Linear_Null::get_U0() noexcept -> arma::vec {
  return U0;
}

auto SKATR_Linear_Null::get_s2() noexcept -> double {
  return s2;
}

auto SKATR_Linear_Null::get_Ux() noexcept -> const arma::mat & {
  return Ux;
}

auto SKATR_Linear_Null::get_coef() noexcept -> const arma::rowvec & {
  return coef;
}

auto SKATR_Linear_Null::get_Y() noexcept -> arma::vec {
  return Y;
}
auto SKATR_Linear_Null::get_X() noexcept -> const arma::mat & {
  return X;
}
