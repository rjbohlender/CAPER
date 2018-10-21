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

  GLM(arma::mat &X, arma::vec &Y, LinkT &link) {
    try{
      beta_ = irls_svdnewton(X, Y);
    } catch(std::exception &e) {
      std::cerr << "IRLS failed; Using gradient descent." << std::endl;
      beta_ = gradient_descent(X, Y);
    }

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
  std::cerr << "Running IRLS.\n";
  const auto tol = 1e-8;
  const auto max_iter = 25;
  auto update = 0.;
  auto iter = 0;

  // Aliases
  arma::mat A = X.t();
  arma::vec &b = Y;

  arma::uword m = A.n_rows;
  arma::uword n = A.n_cols;

  // SVD
  arma::mat U, V;
  arma::vec S;

  arma::svd_econ(U, S, V, A);

  // Matrices and Vectors
  arma::vec t(m, arma::fill::randn);
  arma::vec s(n, arma::fill::randn);
  arma::vec s_old;
  arma::vec weights(m, arma::fill::ones);

  do {
    arma::vec gv = link.linkinv(t);
    arma::vec varg = link.variance(gv);
    arma::vec gvp = link.mueta(t);

    arma::vec z(m, arma::fill::zeros);
    arma::vec W(m, arma::fill::zeros);

    z = t + (b - gv) / gvp;
    W = weights % arma::pow(gvp, 2) / varg;


    s_old = s;

    arma::mat C;
    bool success = arma::chol(C, U.t() * (U.each_col() % W));
    if(!success) {
      std::cerr << z.t();
      std::cerr << W.t();
      std::cerr << gv.t();
      std::cerr << varg.t();
      std::cerr << gvp.t();
    }

    s = solve(arma::trimatl(C.t()), U.t() * (W % z));
    s = solve(arma::trimatu(C), s);

    t = U * s;

    update = arma::norm(s - s_old);

    iter++;
  } while(iter < max_iter && update > tol);

  return (V * (arma::diagmat(1. / S) * (U.t() * t)));
}

#endif //PERMUTE_ASSOCIATE_GLM_HPP
