#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include "frprmn.h"
#include "gauss.h"
#include "reconstructor_GVIR_neto.h"

int main (int argc, char **argv) {
  int n = 500, i, j, nx = 128, ny = 128, k;
  float ftol = 1e-10;
  double **im, *pars, pars2[9];
  int **mask;
  int entropia = 0;
  //  MallaVoronoi *m;
  Reconstructor *r;

  printf ("r = %d\n", r);
  r = leerArchivoEntrada (argv[1]);
  printf ("r = %d\n", r);

  printf ("r->nVis = %d\n\n", r->nVis);
  printf ("ok, %d\n", r->nVis);
  //reinicializar (r, "MEM_NR_002.dat");
  run (r, ftol);

  return 0;
}

