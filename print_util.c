#include "print_util.h"

void printscientificnotation(double var, FILE *fileid, bool printfraction) {
  double varlog10, varfraction;
  varlog10 = floor(log10(var));
  fprintf(fileid, "$");
  if (printfraction) {
    varfraction = var / pow(10.,varlog10);
    fprintf(fileid, "%.2f \\times ", varfraction);
  }
  fprintf(fileid,  "10^{%.0f",varlog10);
  if (fabs(varlog10) < 10)
    fprintf(fileid,  "\\hphantom{0}");
  fprintf(fileid, "}$");
}
