#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include "frprmn.h"
#include "gauss.h"
#include "reconstructor.h"

int main (int argc, char **argv) {
  int n = 3, i, j, nx = 128, ny = 128;
  float ftol = 1e-10;
  double **im, pars[9], pars2[9];
  int **mask;
  MallaVoronoi *m;
  Reconstructor *r;

  r = newReconstructor (argv[1], argv + 2, argc - 2, n);
  //reinicializar (r, "MEM_NR_065.dat");
  //run (r, 0.0, 0, ftol);

  //pars = (double *) malloc (3 * n * sizeof (double));
  for (i = 0; i < 3 * n; i++) {
    //pars[i] = (double) rand() / RAND_MAX;
    pars[i] = i / (3.0 * n);
  }
  printf ("double  = %d\n", sizeof (double));
  printf ("pars    = %d\n", sizeof (pars));
  printf ("double* = %d\n", sizeof (double *));
  printf ("pars2   = %d\n", sizeof (pars2));
  printf ("*pars   = %d\n", sizeof (*pars));
  printf ("        = %d\n", 3 * n * sizeof (double));
  
  srand (0);
  for (i = 0; i < 1e3; i++) {
    //int j1 = 3 * n * ((double) rand() / RAND_MAX);
    int j1 = 3;
    int j2 = 3 * n * ((double) rand() / RAND_MAX);
    //int j2 = 6;
    double aux = pars[j1];
    
    printf ("cambiando %d por %d\n", j1, j2);
    pars [j1] = pars[j2];
    pars [j2] = aux;
    
    printf ("L%d: %g\n", i, L (r->fL, pars, n));
  }

  return 0;
}
