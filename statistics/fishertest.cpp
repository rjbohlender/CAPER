//
// Created by Bohlender,Ryan James on 2019-06-06.
//

#include "fishertest.hpp"
#include <algorithm>
#include <boost/math/distributions/hypergeometric.hpp>

FisherTest::FisherTest(Gene &gene, arma::vec &Y, const std::string &ts)
	: or_(1), p_(1) {
  // Build 2x2 table
  arma::sp_mat X = gene.get_matrix(ts);
  arma::vec collapse(arma::sum(X, 1));
  collapse(arma::find(collapse > 0)).ones();

  double ncase = 2 * arma::accu(Y);
  double ncont = 2 * arma::accu(1 - Y);

#if 0
  case_alt = arma::accu(arma::vec(arma::sum(X, 1) % Y) > 0);
  cont_alt = arma::accu(arma::vec(arma::sum(X, 1) %  (1 - Y)) > 0);
  case_ref = ncase - case_alt;
  cont_ref = ncont - cont_alt;
#else
  case_alt = arma::dot(collapse, Y);
  cont_alt = arma::dot(collapse, (1 - Y));
  case_ref = ncase - case_alt;
  cont_ref = ncont - cont_alt;
#endif
  // Following R's implementation
  // 2x2 table
  //      alt, ref
  //case   x         k
  //cont
  //      m,  n

  // Marginals
  double m = case_alt + cont_alt;
  double n = case_ref + cont_ref;
  double k = case_alt + case_ref;
  double x = case_alt;

  double lo = std::max(0., k - n);
  double hi = std::min(k, m);

  arma::vec support = arma::regspace(lo, hi);
  arma::vec logdc(support.n_elem);

  boost::math::hypergeometric hyper(m, k, m + n);
  arma::uword i = 0;
  for (const auto &v : support) {
	logdc(i) = std::log(boost::math::pdf(hyper, v));
	i++;
  }

  auto dnhyper = [&](double ncp) -> arma::vec {
	arma::vec d = logdc + std::log(ncp) * support;
	d = arma::exp(d - arma::max(d));
	return d / arma::sum(d);
  };

  auto mnhyper = [&](double ncp) -> double {
	if (ncp == 0) {
	  return lo;
	} else if (!std::isfinite(ncp)) {
	  return hi;
	}
	return arma::sum(support % dnhyper(ncp));
  };

  auto pnhyper = [&](double q, double ncp = 1, bool upper_tail = false) -> double {
	if (ncp == 1) {
	  if (upper_tail) {
	    if(x == 0) {
	      return 1;
	    }
		return (boost::math::cdf(boost::math::complement(hyper, x - 1)));
	  } else {
		return (boost::math::cdf(hyper, x));
	  }
	} else if (ncp == 0) {
	  if (upper_tail) {
		return static_cast<double>(q <= lo);
	  } else {
		return static_cast<double>(q >= lo);
	  }
	} else if (!std::isfinite(ncp)) {
	  if (upper_tail) {
		return static_cast<double>(q <= hi);
	  } else {
		return static_cast<double>(q >= hi);
	  }
	}
	if (upper_tail) {
	  return (arma::sum(dnhyper(ncp)(arma::find(support >= q))));
	} else {
	  return (arma::sum(dnhyper(ncp)(arma::find(support <= q))));
	}
  };

  // Calculate P-value Two-Sided
  // if(or_ == 0) {
  // 	p_ = static_cast<double>(x == lo);
  // } else if(!std::isfinite(or_)) {
  // 	p_ = static_cast<double>(x == hi);
  // } else {
  //   double relErr = 1 + 1e-7;
  //   arma::vec d = dnhyper(or_);
  //   p_ = arma::sum(d(arma::find(d <= d(x - lo) * relErr)));
  // }
  p_ = pnhyper(x, or_, true);

  // MLE of the odds ratio
  auto mle = [&](double x) -> double {
    if(x == lo) {
      return 0;
    } else if(x == hi) {
      return arma::datum::inf;
    }
    double mu = mnhyper(1);
	int digits = std::numeric_limits<double>::digits - 3;
	boost::math::tools::eps_tolerance<double> tol(digits);
	boost::uintmax_t max_iter = 1000;
	std::pair<double, double> tmp;
    if(mu > x) {
      tmp = boost::math::tools::toms748_solve([&](double t){ return mnhyper(t) - x; }, 0., 1., tol, max_iter);
	  return (tmp.second + tmp.first) / 2;
    } else if(mu < x) {
	  tmp = boost::math::tools::toms748_solve([&](double t){ return mnhyper(1. / t) - x; }, std::numeric_limits<double>::epsilon(), 1., tol, max_iter);
	  return 1. / ((tmp.second + tmp.first) / 2);
    } else {
      return 1;
    }
  };

  or_ = mle(x);
}

auto FisherTest::get_or() -> double {
  return or_;
}

auto FisherTest::get_pval() -> double {
  return p_;
}

