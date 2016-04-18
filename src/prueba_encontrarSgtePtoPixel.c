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
  int n = 500, i, j, nx = 128, ny = 128, k;
  float ftol = 1e-10;
  double **im, *pars, pars2[9];
  int **mask;
  //  MallaVoronoi *m;
  Reconstructor *r;

  //r = newReconstructor (argv[1], argv + 2, argc - 2, n, 1.0, 0);
  //reinicializar (r, "MEM_NR_305.dat");
  //run (r, ftol);

  PuntoVoronoi *p1 = newPuntoVoronoi(0, 0);
  PuntoVoronoi *p2 = newPuntoVoronoi(0, 0);
  PuntoVoronoi *q1 = newPuntoVoronoi(0, 0);
  PuntoVoronoi *q2 = newPuntoVoronoi(0, 0);
  PuntoVoronoi *pIni = newPuntoVoronoi(0, 0),
    *pFin = newPuntoVoronoi(0, 0);
  double m, dx, dy;
  int tipoLinea, inter, lineaAnterior;

  r = newReconstructor (argv[1], argv + 2, argc - 2, n, 1.0, 0);
  nx = r->fL->fg_image->size[0];
  ny = r->fL->fg_image->size[1];

  for (k = 0; k < 1e5; k++) {
    q1->x = (double) rand() / RAND_MAX * 2 - 0.5;
    q1->y = (double) rand() / RAND_MAX * 2 - 0.5;
    q2->x = (double) rand() / RAND_MAX * 2 - 0.5;
    q2->y = (double) rand() / RAND_MAX * 2 - 0.5;

    printf ("k = %d, q1 = (%g, %g), q2 = (%g, %g)\n", 
	    k, q1->x, q1->y, q2->x, q2->y);

    lineaAnterior = 0;
    m = (q1->y - q2->y) / (q1->x - q2->x);
    if (m > 1) {
      tipoLinea = 3;
    }
    else if (m > 0) {
      tipoLinea = 1;
    }
    else if (m > -1) {
      tipoLinea = 2;
    }
    else {
      tipoLinea = 4;
    }

    if (q1->x < q2->x) {
	if (!(inter = interseccionCuadrado (q1, q2, p1, p2, tipoLinea))) {
	  continue;
	}
    }
    else {
      if (!(inter = interseccionCuadrado (q2, q1, p2, p1, tipoLinea))) {
	continue;
      }
    }
        
    printf ("        p1 = (%g, %g), p2 = (%g, %g) en pixeles\n", 
	    k, p1->x * (nx - 1), p1->y * (ny - 1),
	    p2->x * (nx - 1), p2->y * (ny - 1));

    for (pIni = newPuntoVoronoi (p1->x, p1->y);
	 (pFin->x != p2->x) && (pFin->y != p2->y);) {
      lineaAnterior = encontrarSgtePtoPixel (r->fL, pIni, p2, pFin, lineaAnterior);
      /*
      fprintf (stderr, "      pIni = (%g, %g) pFin =  (%g, %g) (en pixeles)\n", 
	       pIni->x * (nx - 1), pIni->y * (nx - 1), 
	       pFin->x * (ny - 1), pFin->y * (ny - 1));
      */
      i = round (0.5 * (pIni->x + pFin->x) * (nx - 1.0));
      j = round (0.5 * (pIni->y + pFin->y) * (ny - 1.0));    
      
      if (pIni->x * (nx - 1) - (i + 0.5) > 1e-6 ||
	  pIni->x * (nx - 1) - (i - 0.5) < -1e-6 ||
	  pIni->y * (ny - 1) - (j + 0.5) > 1e-6 ||
	  pIni->y * (ny - 1) - (j - 0.5) < -1e-6) {
	fprintf (stderr, "ERROR: pIni = (%g, %g) cae fuera de (%d, %d)\n",
		 pIni->x, pIni->y, i, j);
	exit (1);
      }
      if (pFin->x * (nx - 1) - (i + 0.5) > 1e-6 ||
	  pFin->x * (nx - 1) - (i - 0.5) < -1e-6 ||
	  pFin->y * (ny - 1) - (j + 0.5) > 1e-6 ||
	  pFin->y * (ny - 1) - (j - 0.5) < -1e-6) {
	fprintf (stderr, "ERROR: pFin = (%g, %g) cae fuera de (%d, %d)\n",
		 pFin->x, pFin->y, i, j);
	exit (1);
      }

      if (fabs (m * (pFin->x - p1->x) + p1->y - pFin->y) > 1e-6 ||
	  fabs (m * (pIni->x - p1->x) + p1->y - pIni->y) > 1e-6) {
	fprintf (stderr, "ERROR: punto no pertenece a linea\n");
	fprintf (stderr, "       pIni = (%g, %g), pFin = (%g, %g)\n", 
		 pIni->x, pIni->y, pFin->x, pFin->y);
	fprintf (stderr, "       q1 = (%g, %g), q2 = (%g, %g)\n", 
		 q1->x, q1->y, q2->x, q2->y);
	exit (1);
      }
      if (SQR(pIni->x - pFin->x) + SQR(pIni->y - pFin->y) > 
	  SQR(1.0/(nx - 1.0)) + SQR(1.0/(ny - 1.0))) {
	fprintf (stderr, "ERROR integraldVdx: d(pIni, pFin) muy grande.\n");
	fprintf (stderr, "      pIni = (%g, %g), pFin = (%g, %g)\n", 
		 pIni->x, pIni->y, pFin->x, pFin->y);
	fprintf (stderr, "      d = %.10g vs %.10g\n", 
		 SQR(pIni->x - pFin->x) + SQR(pIni->y - pFin->y), 
		 SQR(1.0/(nx - 1.0)) + SQR(1.0/(ny - 1.0)));
	fprintf (stderr, "      dx = %.10g, dy =  %.10g (en pixeles)\n", 
		 (pFin->x - pIni->x) * (nx - 1), (pFin->y - pIni->y) * (ny - 1));
	fprintf (stderr, "      pIni = (%g, %g) pFin =  (%g, %g) (en pixeles)\n", 
		 pIni->x * (nx - 1), pIni->y * (nx - 1), 
		 pFin->x * (ny - 1), pFin->y * (ny - 1));
	exit (1);    
      }
      
    
      pIni->x = pFin->x;
      pIni->y = pFin->y;
    }

    if (inter && (fabs (m * (q1->x - p1->x) + p1->y - q1->y) > 1e-6 ||
	fabs (m * (q2->x - p1->x) + p1->y - q2->y) > 1e-6)) {
      fprintf (stderr, "ERROR: punto no pertenece a linea\n");
      fprintf (stderr, "       p1 = (%g, %g), p2 = (%g, %g)\n", 
	       p1->x, p1->y, p2->x, p2->y);
      fprintf (stderr, "       q1 = (%g, %g), q2 = (%g, %g)\n", 
	       q1->x, q1->y, q2->x, q2->y);
    }
    
    if (inter && (p1->x > 1 || p1->x < 0 || p2->x > 1 || p2->x < 0 ||
		  p1->y > 1 || p1->y < 0 || p2->y > 1 || p2->y < 0)) {
      fprintf (stderr, "ERROR: punto fuera del cuadrado\n");
      fprintf (stderr, "       p1 = (%g, %g), p2 = (%g, %g)\n", 
	       p1->x, p1->y, p2->x, p2->y);
      fprintf (stderr, "       q1 = (%g, %g), q2 = (%g, %g)\n", 
	       q1->x, q1->y, q2->x, q2->y);
    }
    if (!inter) {
      fprintf (stderr, "NO\n");
    }
    free (pIni);
  }

  return 0;
}
