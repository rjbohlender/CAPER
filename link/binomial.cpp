//
// Created by Bohlender,Ryan James on 10/15/18.
//

#include <algorithm>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>

#include "binomial.hpp"

const std::vector<std::string> Binomial::links {"logit", "probit", "cloglog", "Cauchit", "Log"};

Binomial::Binomial(const std::string &link)
: linkid(check_linkid(link)) {}

arma::vec Binomial::link(arma::mat &X, arma::vec &beta) noexcept {
  switch(linkid) {
	case Binomial::LinkID::Logit:
	  return arma::log(X * beta / (1. - X * beta));
	case Binomial::LinkID::Probit: {
	  arma::vec ret(beta.n_elem);
	  auto vit = ret.begin();
	  boost::math::normal dist(0.0, 1.0);
	  for (auto &v : arma::vec(X * beta)) {
		*vit = boost::math::quantile(dist, v);
		vit++;
	  }
	  return ret;
	}
	case Binomial::LinkID::cloglog:
	  return arma::log(-arma::log((1. - X * beta)));
	case Binomial::LinkID::Cauchit: {
	  arma::vec ret(beta.n_elem);
	  auto vit = ret.begin();
	  boost::math::cauchy dist(0.0, 1.0);
	  for (auto &v : arma::vec(X * beta)) {
		*vit = boost::math::quantile(dist, v);
		vit++;
	  }
	  return ret;
	}
	case Binomial::LinkID::Log:
	  return arma::log(X * beta);
  }
}

arma::vec Binomial::linkinv(arma::mat &X, arma::vec &beta) noexcept {
  switch(linkid) {
  case Binomial::LinkID::Logit:
	return 1. / (1. + arma::exp(-X * beta));
  case Binomial::LinkID::Probit:
    return arma::normcdf(X * beta);
  case Binomial::LinkID::cloglog: {
    arma::vec ret(beta.n_elem);
	auto vit = ret.begin();
	for (auto &v : arma::vec(X * beta)) {
	  *vit = -std::expm1(-std::exp(v));
	  vit++;
	}
	ret = arma::clamp(ret, std::numeric_limits<double>::min(), 1. - std::numeric_limits<double>::min());
	return ret;
  }
  case Binomial::LinkID::Cauchit: {
	arma::vec ret(beta.n_elem);
	auto vit = ret.begin();
	boost::math::cauchy dist(0.0, 1.0);
	for (auto &v : arma::vec(X * beta)) {
	  *vit = boost::math::cdf(dist, v);
	  vit++;
	}
	return ret;
  }
  case Binomial::LinkID::Log:
    return arma::exp(X * beta);
  }
}

arma::vec Binomial::variance(arma::vec &mu) noexcept {
  return mu % (1. - mu);
}

arma::vec Binomial::mueta(arma::vec &eta) noexcept {
  switch(linkid) {
  case LinkID::Logit:
	return arma::exp(eta) / arma::pow((arma::exp(eta) + 1.), 2);
  case LinkID::Probit:
    return arma::normpdf(eta);
  case LinkID::cloglog: {
    eta = arma::clamp(eta, arma::min(eta), 700);
    arma::vec ret = arma::exp(eta) * arma::exp(-arma::exp(eta));
	return arma::clamp(ret, std::numeric_limits<double>::min(), arma::max(ret));
  }
  case LinkID::Cauchit: {
    boost::math::cauchy dist(0., 1.);
	arma::vec ret(eta.n_elem);
	auto vit = ret.begin();
	for (auto &v : eta) {
	  *vit = boost::math::pdf(dist, v);
	  vit++;
	}
	ret = arma::clamp(ret, std::numeric_limits<double>::min(), 1. - std::numeric_limits<double>::min());
  }
  case LinkID::Log:
    return arma::clamp(arma::exp(eta), std::numeric_limits<double>::min(), 1. - std::numeric_limits<double>::min());
  }
}

Binomial::LinkID Binomial::check_linkid(const std::string &link) {
  auto ok = std::find_if(links.cbegin(), links.cend(), linkid);
  if(ok == links.cend())
    throw(std::logic_error("Wrong link argument to Binomial."));

  if(link == "logit") {
    return Binomial::LinkID::Logit;
  } else if(link == "probit") {
    return Binomial::LinkID::Probit;
  } else if(link == "cloglog") {
    return Binomial::LinkID::cloglog;
  } else if(link == "cauchit") {
    return Binomial::LinkID::Cauchit;
  } else {
    return Binomial::LinkID::Log;
  }
}
