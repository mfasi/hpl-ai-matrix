#include "matgen.h"

/*******************
 * Generate matrix *
 *******************/

// Generate pivoted matrix with pre-assigned inf-norm condition number.
void generateA(matrixtype *A,
               const MKL_INT *descA,
               const MKL_INT nA,
               const double alpha,
               const double beta,
               const MKL_INT mympirank,
               const MKL_INT nprows,
               const MKL_INT npcols) {

  // Parse descriptor.
  MKL_INT ctxt = descA[1]; // BLACS context.
  MKL_INT M = descA[2];    // Rows of global matrix.
  MKL_INT N = descA[3];    // Columns of global matrix.
  MKL_INT Mb = descA[4];   // Blocking factor for rows.
  MKL_INT Nb = descA[5];   // Blocking factor for columns.
  MKL_INT mA = descA[8];
  MKL_INT ONE = 1;

  MKL_INT prow, pcol;
  blacs_pcoord(&ctxt, &mympirank, &prow, &pcol);

  double ab = alpha * beta;
  double a = -alpha;
  double b = -beta;
  MKL_INT iloc, jloc, jdisp, i, j, k1, k2;
  jloc = 0;
  j = pcol*Nb+1;
  while (jloc<nA){
    for (k1=0; k1<Nb && j<=N; k1++) {
      jdisp = jloc*mA;
      i = prow*Mb+1;
      iloc = 0;
      while (iloc<mA) {
        for (k2=0; k2<Mb && i<=M; k2++) {
          if (i > j)
            A[jdisp+iloc] = (matrixtype)(a + (j-1) * ab);
          else if (i == j)
            A[jdisp+iloc] = (matrixtype)(1 + (i-1) * ab);
          else // (j < i)
            A[jdisp+iloc] = (matrixtype)(b + (i-1) * ab);
          iloc++;
          i++;
        }
        i+=Mb*(nprows-1);
      }
      jloc++;
      j++;
    }
    j+=Nb*(npcols-1);
  }
}





/************************************
 * Maximum lower triangular element *
 ************************************/

// Compute max element (in magnitude) and max error of lower triangular factor.
void maxlowertri(matrixtype *A,
                 const double alpha,
                 const MKL_INT *descA,
                 const MKL_INT nA,
                 const MKL_INT mympirank,
                 const MKL_INT nprows,
                 const MKL_INT npcols,
                 double *maxrelerror) {

  // Parse descriptor.
  MKL_INT ctxt = descA[1]; // BLACS context.
  MKL_INT M = descA[2];    // Rows of global matrix.
  MKL_INT N = descA[3];    // Columns of global matrix.
  MKL_INT Mb = descA[4];   // Blocking factor for rows.
  MKL_INT Nb = descA[5];   // Blocking factor for columns.
  MKL_INT mA = descA[8];
  MKL_INT ONE = 1;

  MKL_INT prow, pcol;
  blacs_pcoord(&ctxt, &mympirank, &prow, &pcol);

  MKL_INT iloc, jloc, jdisp, i, j, k1, k2;
  jloc = 0;
  j = pcol*Nb+1;
  *maxrelerror = -INFINITY;
  double candidate;
  while (jloc<nA){
    for (k1=0; k1<Nb && j<=N; k1++) {
      jdisp = jloc*mA;
      i = prow*Mb+1;
      iloc = 0;
      while (iloc<mA) {
        for (k2=0; k2<Mb && i<=M; k2++) {
          if (i > j) {
            candidate = fabs(A[jdisp+iloc] + alpha);
            if( candidate / alpha > *maxrelerror)
              *maxrelerror = (double)candidate;
          }
          iloc++;
          i++;
        }
        i+=Mb*(nprows-1);
      }
      jloc++;
      j++;
    }
    j+=Nb*(npcols-1);
  }

}
