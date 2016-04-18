#ifndef GAUSS_GCV
#define GAUSS_GCV

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif  

typedef struct {
  int n;        // Numero de sumas.
  double media;
  double sigma;
  double gauss_add, gauss_fac;
} Gaussiana;

Gaussiana *newGaussiana (int n, double media, 
			 double sigma, int seed);
double gauss (Gaussiana *g);

#ifdef __cplusplus
}
#endif  

#endif
