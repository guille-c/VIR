#define NRANSI
#include "nrutil.h"
#include "linmin.h"
#include <stdarg.h>
#include <stdio.h>
#define TOL 2.0e-4

int ncom;
float *pcom,*xicom,(*nrfunc)(float []);

void imprimirLogLinmin (char* s, ...) {
  va_list ap;
  FILE* archivoLog = fopen ("reconstructor.log", "a");

  va_start (ap, s);
  vfprintf (archivoLog, s, ap);
  va_end (ap);

  fclose (archivoLog);
}

void linmin(float p[], float xi[], int n, float *fret, float (*func)(float []))
{
	float brent(float ax, float bx, float cx,
		float (*f)(float), float tol, float *xmin);
	float f1dim(float x);
	void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
		float *fc, float (*func)(float));
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;
	float promedio = 0; // GUILLE 25/04/2006

	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
		promedio += fabs (xi[j]); // GUILLE 25/04/2006
		//printf ("linmin: xi[%d] = %g, promedio = %g\n", j, xi[j], promedio);
	}
	promedio /= n; // GUILLE 25/04/2006
	ax=0.0;
	xx=1.0; // punt
	xx = 1e-3; // GUILLE 01/06/2007
	//xx = 1e4; // helix CBI
	//xx = 1e-1; //gaussianas 4
	//xx = 1e-5; //gaussianas 1
	//xx = 1e-5; // L1622
	//xx = 1e-6;
	//xx = 1e-8; // ROPHW
	//xx = 1e-17;
	//xx = 0.1 / promedio; // GUILLE 25/04/2006
	printf ("linmin: xx = %g, promedio = %g, n = %d\n", xx, promedio, n);
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	printf ("-------finalizado mnbrak a = %g, b = %g, c = %g\n", ax, xx, bx);
	imprimirLogLinmin ("-------finalizado mnbrak a = %g, b = %g, c = %g\n", ax, xx, bx);
	imprimirLogLinmin ("                        fa = %g, fb = %g, fc = %g\n", fa, fx, fb);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	imprimirLogLinmin ("-------finalizado brent xmin = %g, fret = %g\n",
		     xmin, *fret);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
		//imprimirLogLinmin ("       p[%d] = %g\n", j, xi[j]);
	}
	imprimirLogLinmin ("\n");
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}
#undef TOL
#undef NRANSI
