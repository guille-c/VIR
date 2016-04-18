#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include <gsl/gsl_multimin.h>

funcL *f;


double
my_f (const gsl_vector *v, void *params)
{
  double ret = 0;
  int n = ((int *)params)[0], i;
  double *pars = malloc (3 * n * sizeof(double));
  
  for (i = 0; i < n; i++) {
    pars[3 * i] = gsl_vector_get(v, 3 * i);
    pars[3 * i + 1] = gsl_vector_get(v, 3 * i + 1);
    pars[3 * i + 2] = gsl_vector_get(v, 3 * i + 2);
  }
  ret = L(f, pars, n);
  free (pars);
  return ret;
}

int main(int argc, char **argv) {
  int n = 500, par[1] = {n};
  size_t iter = 0, np = n * 3, i;
  int status;
  double ret;
  gsl_vector *x;

  f = newFuncL(argv[1], argv + 2, argc - 2);
  x = gsl_vector_alloc (np);
  for (i = 0; i < n; i++) {
    gsl_vector_set(x, 3 * i, (double) random() / RAND_MAX);
    gsl_vector_set(x, 3 * i + 1, (double) random() / RAND_MAX);
    gsl_vector_set(x, 3 * i + 2, 0.0);
  }

  for (i = 0; i < 1e20; i++) {
    ret = my_f(x, (void *)par);
  }
  return 0;
}
