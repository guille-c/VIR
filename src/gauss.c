#include "gauss.h"

Gaussiana *newGaussiana(int n, double media,
			double sigma, int seed) {
  Gaussiana *g = (Gaussiana *) malloc (sizeof (Gaussiana));
  g->n = n;
  g->media = media;
  g->sigma = sigma;
  g->gauss_add = sqrt(3 * n);
  g->gauss_fac = 2 * g->gauss_add / ( (double )n * RAND_MAX);
  srand (seed);
  return g;
}

double gauss(Gaussiana *g)
{
  double sum;
  int i;
  
  sum = 0;
  for (i = 0; i < g->n; i++)
    sum += rand();
  return (g->gauss_fac * sum - g->gauss_add) * g->sigma 
    + g->media;
}
