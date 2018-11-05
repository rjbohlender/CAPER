//
// Created by Bohlender,Ryan James on 10/16/18.
//

#ifndef PERMUTE_ASSOCIATE_GLM_HPP
#define PERMUTE_ASSOCIATE_GLM_HPP

#include <armadillo>

#include "../link/family.hpp"

template <typename LinkT>
struct GLM {
  LinkT link;
  arma::vec beta_; // coefficients
  arma::vec mu_;   // fitted.values
  arma::vec eta_;  // linear.predictors
  double dev_;

  GLM(arma::mat &X, arma::vec &Y, LinkT &link) {
    //try{
      beta_ = irls_svdnewton(X, Y);
    //} catch(std::exception &e) {
      //std::cerr << "IRLS failed; Using gradient descent." << std::endl;
      //beta_ = gradient_descent(X, Y);
      //std::cerr << e.what();
    //}

    arma::mat A = X.t();
    mu_ = link.linkinv(A, beta_);
    eta_ = link.link(A, beta_);
  }

  // Algorithms for finding the optimum
  auto gradient_descent(arma::mat &X, arma::colvec &Y) -> arma::vec;
  auto irls_svdnewton(arma::mat &X, arma::colvec &Y) -> arma::vec;
};

template<typename LinkT>
auto GLM<LinkT>::gradient_descent(arma::mat &X, arma::colvec &Y) -> arma::vec {
  std::cerr << "Running gradient descent.\n";
  auto iterations = 0ull;
  const auto max_iter = 1000000ull;
  auto alpha = 0.005; // Learning rate
  auto tol = 1e-10;
  auto m = static_cast<double>(X.n_cols);
  arma::mat A = X.t();
  auto b = arma::vec(A.n_cols, arma::fill::randn);
  auto grad = b;
  do {
    // Vectorized update
    grad = alpha * (A.t() * (link.linkinv(A, b) - Y)) / m;
    b -= grad;


    iterations++;
  } while(iterations < max_iter && arma::norm(grad) > tol);
  return b;
}

template<typename LinkT>
auto GLM<LinkT>::irls_svdnewton(arma::mat &X, arma::colvec &Y) -> arma::vec {
  arma::wall_clock timer;
  timer.tic();
  const auto tol = 1e-8;
  const auto max_iter = 25;
  auto iter = 0;

  arma::uword m = X.n_cols;
  arma::uword n = X.n_rows;

  // SVD
  arma::mat U, V;
  arma::vec S;

  arma::svd_econ(U, S, V, X.t(), "both", "std");

  // Matrices and Vectors
  arma::vec eta(m, arma::fill::randn);
  arma::vec s(n, arma::fill::randn);
  arma::vec weights(m, arma::fill::ones);

  dev_ = 0;
  double devold;

  do {
    devold = dev_;

    arma::vec g = link.linkinv(eta);
    arma::vec varg = link.variance(g);
    arma::vec gprime = link.mueta(eta);

    arma::vec z(m, arma::fill::zeros);
    arma::vec W(m, arma::fill::zeros);

    z = eta + (Y - g) / gprime;
    W = weights % arma::pow(gprime, 2) / varg;

    arma::mat C;
    arma::mat UWU = U.t() * (U.each_col() % W);
    bool success = arma::chol(C, UWU);
    arma::vec scale = UWU.diag();
    arma::uword tries = 0;
    while(!success && tries < 10) {
      UWU.diag() += arma::mean(scale) + 1e-5;
      success = arma::chol(C, UWU);
      tries++;
    }

    s = solve(arma::trimatl(C.t()), U.t() * (W % z));
    s = solve(arma::trimatu(C), s);

    eta = U * s;

    iter++;

    dev_ = arma::sum(link.dev_resids(Y, g, weights));

  } while(iter < max_iter && std::abs(dev_ - devold) / (0.1 + std::abs(dev_)) > tol);

  return (V * (arma::diagmat(1. / S) * (U.t() * eta)));
}

#endif //PERMUTE_ASSOCIATE_GLM_HPP
