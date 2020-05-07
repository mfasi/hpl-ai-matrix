#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include <stdbool.h>
#include <math.h>

#include "matnorm.h"
#include "print_util.h"

/*************
 * RUN TESTS *
 *************/
int main(int argc, char **argv) {

  // Parse input arguments.
  static int debugint = false;
  static struct option
    long_opts[] = {{"verbose", no_argument, &debugint, 1},
                   {"nprows", required_argument, 0, 'm'},
                   {"npcols", required_argument, 0, 'n'},
                   {0, 0, 0, 0}
  };
  int option_ind = 0;
  int opt;
  while ((opt = getopt_long(argc, argv, "v", long_opts, &option_ind)) != -1)
    {
      switch (opt) {
      case 0:
        break;
      case 'v':
        debugint = true;
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

  // Declarations.
  double alpha, beta, alphabeta, fres;
  double iopt, jopt;

  char outfilenamebeta [] = "./results-beta.dat";
  char outfilenameres [] = "./results-fres.dat";
  char outfilenamemin [] = "./results-min.dat";
  char outfilenamemax [] = "./results-max.dat";
  FILE *outfilebeta = fopen(outfilenamebeta, "w");
  FILE *outfileres = fopen(outfilenameres, "w");
  FILE *outfilemin = fopen(outfilenamemin, "w");
  FILE *outfilemax = fopen(outfilenamemax, "w");
  assert(outfilebeta != NULL);
  assert(outfilemin != NULL);
  assert(outfilemax != NULL);

  size_t i, j, n;
  size_t sizes [] = {1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10};
  size_t nsizes = 9;
  double ratio = 0.5, largest_elem, smallest_elem;
  double kappa;
  double kappas [] = {1e2, 1e4, 1e6, 1e8, 1e10};
  size_t nkappas = 5;

  for (i=0; i<nsizes; i++) {
    n = sizes[i];

    printscientificnotation(n, outfilebeta, false);
    printscientificnotation(n, outfileres, false);
    printscientificnotation(n, outfilemin, false);
    printscientificnotation(n, outfilemax, false);

    for (j=0; j<nkappas; j++) {
      kappa = kappas[j];

      // Find right alpha and beta.
      beta = findparameters(n, ratio, kappa);
      fres = fabs(condA(n,ratio*beta,beta)-kappa)/kappa;
      alpha = ratio * beta;
      alphabeta = alpha * beta;

      // Largest element of generated matrix.
      largest_elem = 1 + (n-1) * alphabeta;

      // Smallest element of generated matrix.
      jopt = fmin(round(1/beta+1), n);
      iopt = fmin(round(1/alpha+1), n);
      smallest_elem = fmin(fabs(-alpha + (jopt-1)*alpha*beta),
                           fabs(-beta + (iopt-1)*alpha*beta));

      fprintf(outfilebeta, " & ");
      printscientificnotation(beta, outfilebeta, true);
      fprintf(outfileres, " & ");
      printscientificnotation(fres, outfileres, true);
      fprintf(outfilemin, " & ");
      printscientificnotation(smallest_elem, outfilemin, true);
      fprintf(outfilemax, " & ");
      printscientificnotation(largest_elem-1, outfilemax, true);

      printf("%.2e ", fres);
      printf(" & ");

    }
    fprintf(outfilebeta, "\\\\\n");
    fprintf(outfileres, "\\\\\n");
    fprintf(outfilemin, "\\\\\n");
    fprintf(outfilemax, "\\\\\n");

    printf("\\\\\n");
  }

  // Free resources.
  fclose(outfilebeta);
  fclose(outfileres);
  fclose(outfilemin);
  fclose(outfilemax);

  return 0;
}
