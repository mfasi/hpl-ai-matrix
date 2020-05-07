#ifndef GENMAT_H_
#define GENMAT_H_

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <mpi.h>
#include <mkl.h>
#include <mkl_blacs.h>
#include <mkl_pblas.h>
#include <mkl_scalapack.h>

typedef float matrixtype;

void generateA(matrixtype *, const MKL_INT *, const MKL_INT,
               const double, const double,
               const MKL_INT, const MKL_INT, const MKL_INT);
void maxlowertri(matrixtype *, const double, const MKL_INT *, const MKL_INT,
                 const MKL_INT, const MKL_INT, const MKL_INT,
                 double *);

#endif // GENMAT_H_
