#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include "frprmn.h"
//#include "mcheck.h"

int main (int argc, char **argv) {
  int n = 100, i, indice = 9, t;
  funcL *fL;
  FILE *out;
  double *pars, *grad, *grad_aprox, delta;

  //mtrace();
  
  if (argc < 3) {
    printf ("ERROR usage: comprobacion_gradiente archivo.fits archivo.sub malla.dat\n");
    exit (1);
  }
  
  out = fopen ("gradiente_aprox_100.dat", "w");
  pars = (double *) malloc (3 * n * sizeof(double));
  for (i = 0; i < n; i++) {
    pars[3 * i] = (double) random() / RAND_MAX;
    pars[3 * i + 1] = (double) random() / RAND_MAX;
    pars[3 * i + 2] = (double) random() / RAND_MAX;
  }

  grad       = (double *) malloc (3 * n * sizeof(double));
  grad_aprox = (double *) malloc (3 * n * sizeof(double));
  fL = newFuncL (argv[1], argv + 2, argc - 3, n);
  
  t = time(0);
  //for (i = 0; i < 10; i++) {
    dL (fL, pars, grad, n, 0);
    //}
  printf ("tiempo exacto: %d\n", time (0) - t);
  t = time(0);
  //for (i = 0; i < 10; i++) {
    dL (fL, pars, grad_aprox, n, 1);
    //}
  printf ("tiempo aprox: %d\n", time (0) - t);

  for (i = 0; i < n; i++) {
    fprintf (out, "%.15g\t%.15g\n", grad[i], grad_aprox[i]);
  }
  fclose (out);  
  
  eliminarFuncL (fL);
  free (pars);
  free (grad);
  
  return 0;
}


