//
// Created by Bohlender,Ryan James on 8/22/18.
//

#ifndef PERMUTE_ASSOCIATE_QRFC2_HPP
#define PERMUTE_ASSOCIATE_QRFC2_HPP

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif

#include <armadillo>
#include <csetjmp>

class QFC {
public:
  QFC(std::vector<double> &lb1, std::vector<double> &nc1, std::vector<int> &n1, double sigma, double c1, int lim1, double acc);

  double get_res();
  int get_fault();

private:
  // Field definitions
  double sigsq;
  double lmax;
  double lmin;
  double mean;
  double c;
  double intl;
  double ersm;
  double acc;
  double res;

  int count;
  int r;
  int lim;
  int ifault;

  bool ndtsrt;
  bool fail;

  std::vector<int> n;
  std::vector<int> th;

  std::vector<double> lb;
  std::vector<double> nc;
  std::vector<double> trace;

  jmp_buf env;

  static double exp1(double x);
  static double square(double x);
  static double cube(double x);
  static double log1(double x, bool first);

  void order();
  double errbd(double u, double *cx);
  double ctff(double accx, double *upn);
  double truncation(double u, double tausq);
  void integrate(int nterm, double interv, double tausq, bool mainx);
  void findu(double *utx, double accx);
  double cfe(double x);

  void counter();


  // Static const variables
  static constexpr double pi = 3.14159265358979;
  static constexpr double log28 = 0.0866; // log(2.0) / 8.0

};

#endif //PERMUTE_ASSOCIATE_QRFC2_HPP
