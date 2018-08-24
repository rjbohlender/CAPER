//
// Created by Bohlender,Ryan James on 8/17/18.
//

#ifndef PERMUTE_ASSOCIATE_SKAT_HPP
#define PERMUTE_ASSOCIATE_SKAT_HPP

extern "C" {

void qfc_1(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res);
void qfc(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res);

};


#endif //PERMUTE_ASSOCIATE_SKAT_HPP
