#ifndef VIRPP_FRPRMN
#define VIRPP_FRPRMN

#include <math.h>
#include "nrutil.h"

#ifdef __cplusplus
extern "C" {
#endif  

#define NRANSI
#define ITMAX 200
#define EPS 1.0e-10
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	    float (*func)(float[]), void (*dfunc)(float [], float []));

#ifdef __cplusplus
}
#endif

#endif
