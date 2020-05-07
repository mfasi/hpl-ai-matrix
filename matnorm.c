#include "matnorm.h"

/*************************************
 * Compute norm and condition number *
 *************************************/

// Compute infinity norm using explicit formula.
double normA(size_t n, double a, double b) {
  double tmp = floor((1.+b)/b);
  double  k = tmp<n-1?tmp:n-1;
  double iprime = floor(1./a);
  double gamma1 = 1. + (n-1) * b;
  double max12;
  if (iprime < n) {
    double gamma2 = 1. + (2.*k-iprime+1.)*a +
      (n-iprime)*b*(-k*k + k + 3.*iprime*(iprime-1.)/2. - n*iprime + n)*a*b;
    max12 = gamma1 > gamma2?gamma1:gamma2;
  } else {
    max12 = gamma1;
  }
  double gamma3 = 1. + (2.*k-n+1)*a + (-k*k + k + n*(n-1.)/2.)*a*b;
  return max12 > gamma3?max12:gamma3;
}

// Compute infinity norm of the inverse using explicit formula.
double normAinv(size_t n, double a, double b) {
  double r = (1.+a)*(1.+b);
  double delta1 = 1.+(1.+a)*b*(1.-pow(r,n-1))/(1.-r);
  double delta2 = pow(1.+a,n-1.);
  return delta1>delta2?delta1:delta2;
}

// Compute infinity-norm condition number using explicit formulae.
double condA(size_t n, double a, double b) {
  return normA(n, a, b) * normAinv(n, a, b);
}





/*******************/
/* Find parameters */
/*******************/

// Function used for the optimization.
double optfun(double b, void *params) {
  struct optfun_params *p;
  p = (struct optfun_params *)params;
  return condA(p->n, p->ratio * b, b) - p->kappa;
}

double findparameters(size_t n, double ratio, double kappa) {
  // Initialize GSL function.
  gsl_function gslfun;
  struct optfun_params params = {n, ratio, kappa};
  gslfun.function = &optfun;
  gslfun.params = &params;

  // Set limits so that 0 < a < 1
  double b_lo = DBL_EPSILON; // *n/1000;
  /* double b_hi = 10*sqrt(sqrt(sqrt(sqrt(kappa)))) / n; */
  double b_hi = 1 / ratio;
  while (!isfinite(condA(n,ratio*b_hi,b_hi)))
    b_hi /= 2;

  // Allocate an instance of Brent's solver.
  gsl_root_fsolver *solver;
  solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  gsl_root_fsolver_set(solver, &gslfun, b_lo, b_hi);

  const size_t MAXITER = 100;
  const double TOL = pow(2,-52);
  size_t i;
  int status = GSL_CONTINUE;
  double approx = 0, bound_lo, bound_hi;
  for (i=0; i < MAXITER && status == GSL_CONTINUE; i++) {
    // Perform one step and check convergence.
    status = gsl_root_fsolver_iterate(solver);
    if (status != GSL_SUCCESS) { // error
      break;
    } else {                     // if status == GSL_SUCCESS exits.
      approx = gsl_root_fsolver_root(solver);
      bound_lo = gsl_root_fsolver_x_lower(solver);
      bound_hi = gsl_root_fsolver_x_upper(solver);
      status = gsl_root_test_interval(bound_lo, bound_hi, 0, TOL);
    }
  }
  assert(fabs(condA(n,ratio*approx,approx)-kappa)/kappa < sqrt(FLT_EPSILON)/2);

  // Deallocate the solver.
  gsl_root_fsolver_free(solver);

  // Return approximations, and print warning if needed.
  if (status == GSL_SUCCESS) {
    return approx;
  } else {
    printf("Solver did not converge to tolerance %.2e after %lu iterations.\n",
           TOL, MAXITER);
    return approx;
  }
}
