#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include <stdbool.h>
#include <math.h>

#include <mpi.h>
#include <mkl.h>
#include <mkl_blas.h>
#include <mkl_blacs.h>
#include <mkl_pblas.h>
#include <mkl_scalapack.h>

#include "print_util.h"
#include "matgen.h"
#include "matnorm.h"

/*************
 * RUN TESTS *
 *************/
int main(int argc, char **argv) {

  // Initialize MPI.
  MPI_Init(&argc, &argv);

  // Parse input arguments.
  static int debugint = false;

  // Parse input arguments.
  MKL_INT nprows = 1, npcols = 1;
  static struct option
    long_opts[] = {{"verbose", no_argument, &debugint, 1},
                   {"nprows", required_argument, 0, 'm'},
                   {"npcols", required_argument, 0, 'n'},
                   {0, 0, 0, 0}
  };
  int option_ind = 0;
  int opt;
  while ((opt = getopt_long(argc, argv, "M:m:n:k:", long_opts, &option_ind)) != -1)
    {
      switch (opt) {
      case 0:
        break;
      case 'm':
        nprows = atoi(optarg);
        assert(nprows>0);
        break;
      case 'n':
        npcols = atoi(optarg);
        assert(npcols>0);
        break;
      case '?':
        printf("?\n");
        printf("Unrecognized option %c\n", optopt);
        break;
      case ':':
        printf(":\n");
        printf("Option %c requires an argument\n", optopt);
        break;
      default:
        abort();
      }
    }
  bool debug = (bool)debugint;

  // Initialize BLACS.
  MKL_INT mympirank_blacs, nmpiprocs_blacs;
  int mympirank, nmpiprocs;
  MKL_INT ctxt, prow, pcol;
  MKL_INT info;
  MKL_INT ONE = 1, MINUSONE = -1, ZERO=0;

  blacs_pinfo(&mympirank_blacs, &nmpiprocs_blacs);
  mympirank = (int)mympirank_blacs;
  nmpiprocs = (int)nmpiprocs_blacs;
  blacs_get(&MINUSONE, &ZERO, &ctxt);
  blacs_gridinit(&ctxt, "R", &nprows, &npcols);
  blacs_gridinfo(&ctxt, &nprows, &npcols, &prow, &pcol);


  // Declarations.
  double alpha, beta;
  double tstart, tend, telapsed;
  double tfzero, tmatgen, tlu;

  MKL_INT M;
  MKL_INT mb = 2;
  MKL_INT nb = 2;
  MKL_INT descA[9];


  char outfilenamestability [50];
  char outfilenametiming [50];
  sprintf(outfilenamestability, "./results-stability-%d.dat", nmpiprocs);
  sprintf(outfilenametiming, "./results-timing-%d.dat", nmpiprocs);

  FILE *outfilestability = NULL;
  FILE *outfiletiming = NULL;

  if (mympirank == 0) {
    outfilestability = fopen(outfilenamestability, "w");
    outfiletiming = fopen(outfilenametiming, "w");
    assert(outfilestability != NULL);
    assert(outfiletiming != NULL);
  }

  double ratio = 0.1;

  size_t sizes [] = {1000,2000,5000,10000,20000,50000,100000,200000,};
  double kappa;
  double kappas [] = {1e3, 1e6};

  size_t i, j; //, nsizes = 8, nkappas = 2;
  size_t nsizes = 8, nkappas = 2;

  for (i=0; i<nsizes; i++) {

    M = sizes[i];
    if (mympirank == 0) {
      fprintf(stdout, "%6lld   ", M);
      printscientificnotation(M, outfilestability, true);
      printscientificnotation(M, outfiletiming, true);
    }

    // Initialize test matrix and pivoting array
    MKL_INT mA = numroc(&M, &mb, &prow, &ZERO, &nprows);
    MKL_INT nA = numroc(&M, &nb, &pcol, &ZERO, &npcols);
    descinit(descA, &M, &M, &mb, &nb, &ZERO, &ZERO, &ctxt, &mA, &info);
    matrixtype *A = (matrixtype *)malloc(mA*nA*sizeof(matrixtype));
    MKL_INT *ipiv = (MKL_INT *)malloc((mA+mb)*sizeof(MKL_INT));

    for (j=0; j<nkappas; j++) {

      kappa = kappas[j];

      if (mympirank == 0) {
        tstart = MPI_Wtime();
        beta = findparameters(M, ratio, kappa);
        tend = MPI_Wtime();
        tfzero = tend - tstart; // Defined only on root process.
      }
      MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      alpha = ratio * beta;

      // Generate test matrix.
      tstart = MPI_Wtime();
      generateA(A, descA, nA, alpha, beta, mympirank, nprows, npcols);
      tend = MPI_Wtime();
      telapsed = tend - tstart;
      MPI_Reduce(&telapsed, &tmatgen, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

      // Compute the LU decomposition.
      tstart = MPI_Wtime();
      psgetrf(&M, &M, A, &ONE, &ONE, descA, ipiv, &info);
      tend = MPI_Wtime();
      telapsed = tend - tstart;
      MPI_Reduce(&telapsed, &tlu, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

      // Find value of largest entry in lower triangular factor.
      double local_maxrelerror, global_maxrelerror;
      maxlowertri(A, alpha, descA, nA, mympirank, nprows, npcols,
                  &local_maxrelerror);
      MPI_Reduce(&local_maxrelerror, &global_maxrelerror, 1, MPI_DOUBLE, MPI_MAX, 0,
                 MPI_COMM_WORLD);

      // Write to file. Format:
      // order
      // alpha beta maxrelerr
      // t_fzero t_matgen t_lu
      if (mympirank == 0) {
        fprintf(stdout, "%.2e %.2e    %.1e %.1e %.1e      ",
                beta, global_maxrelerror,
                tfzero, tmatgen, tlu);

        fprintf(outfilestability, " & ");
        printscientificnotation(beta, outfilestability, true);
        fprintf(outfilestability, " & ");
        printscientificnotation(global_maxrelerror, outfilestability, true);
        fflush(outfilestability);

        fprintf(outfiletiming, " & ");
        printscientificnotation(tfzero, outfiletiming, true);
        fprintf(outfiletiming, " & ");
        printscientificnotation(tmatgen, outfiletiming, true);
        fprintf(outfiletiming, " & ");
        printscientificnotation(tlu, outfiletiming, true);
        fflush(outfiletiming);
      }
    }

    // Deallocate memory.
    free(A);
    free(ipiv);

    if (mympirank == 0) {
      fprintf(stdout, "\n");
      fprintf(outfilestability, "\\\\\n");
      fprintf(outfiletiming, "\\\\\n");
    }

  }

  // Free resources.
  blacs_gridexit(&ctxt);
  MPI_Finalize();
  if (mympirank == 0) {
    fclose(outfilestability);
    fclose(outfiletiming);
  }

  return 0;
}
