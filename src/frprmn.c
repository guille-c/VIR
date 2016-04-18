#include "frprmn.h"
#include <stdio.h>

// GUILLE 31/05/2007
void imprimir (float fp, float fr, float *g, int n){
  FILE *archivo = fopen ("frprmin.log", "a");
  int i;

  fprintf (archivo, "fp = %g\n", fp);
  fprintf (archivo, "fr = %g\n", fr);
  fprintf (archivo, "\n");
  for (i = 1; i <= n; i++){
    fprintf (archivo, "grad[%d] =\t%g\n", i, g[i]);
  }
  fprintf (archivo, "\n");
  fclose (archivo);
}

void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []))
{
	void linmin(float p[], float xi[], int n, float *fret,
		float (*func)(float []));
	int j,its;
	float gg,gam,fp,dgg;
	float *g,*h,*xi;

	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		linmin(p,xi,n,fret,func);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
		  (*dfunc)(p,g);
		  imprimir (*fret, fp, g, n);
			FREEALL
			  //printf ("Fin 1 frprmn\n");
			return;
		}
		fp=(*func)(p);
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
		//if (gg <=  ftol) { // GUILLE 31/05/2007
		  imprimir (*fret, fp, g, n);
			FREEALL
			  //printf ("Fin 2 frprmn\n");
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	fprintf(stderr, "Too many iterations in frprmn"); // GUILLE 23/12/2006
	//nrerror("Too many iterations in frprmn"); // GUILLE 23/12/2006
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI
