#ifndef MATNORM_H_
#define MATNORM_H_

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

double normA(size_t, double, double);
double normAinv(size_t, double, double);
double condA(size_t, double, double);

struct optfun_params {
  size_t n;
  double ratio; // a = ratio * b
  double kappa;
};

double optfun(double b, void *params);
double findparameters(size_t , double, double);

#endif // MATNORM_H_
