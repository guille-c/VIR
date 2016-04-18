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
  int n = 500, i, j, nx = 128, ny = 128;
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
  double m, dx, dy;
  int tipoLinea, inter;
  for (i = 0; i < 1e5; i++) {
    p1->x = (double) rand() / RAND_MAX * 2 - 0.5;
    p1->y = (double) rand() / RAND_MAX * 2 - 0.5;
    p2->x = (double) rand() / RAND_MAX * 2 - 0.5;
    p2->y = (double) rand() / RAND_MAX * 2 - 0.5;

    printf ("i = %d, p1 = (%g, %g), p2 = (%g, %g)\n", 
	    i, p1->x, p1->y, p2->x, p2->y);

    m = (p1->y - p2->y) / (p1->x - p2->x);
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
  
    if (p1->x < p2->x) {
      if (inter = interseccionCuadrado (p1, p2, q1, q2, tipoLinea)) {
	if (q1->x < p1->x - 1e-6 || q2->x < p1->x - 1e-6) {
	  fprintf (stderr, "ERROR: punto fuera de la arista\n");
	  fprintf (stderr, "       p1 = (%g, %g), p2 = (%g, %g)\n", 
		   p1->x, p1->y, p2->x, p2->y);
	  fprintf (stderr, "       q1 = (%g, %g), q2 = (%g, %g)\n", 
		   q1->x, q1->y, q2->x, q2->y);
	}
	if (q2->x > p2->x + 1e-6 || q1->x > p2->x + 1e-6) {
	  fprintf (stderr, "ERROR: punto fuera de la arista\n");
	  fprintf (stderr, "       p1 = (%g, %g), p2 = (%g, %g)\n", 
		   p1->x, p1->y, p2->x, p2->y);
	  fprintf (stderr, "       q1 = (%g, %g), q2 = (%g, %g)\n", 
		   q1->x, q1->y, q2->x, q2->y);
	}      
      }
    }
    else {
      if (inter = interseccionCuadrado (p2, p1, q1, q2, tipoLinea)) {
	if (q2->x > p1->x + 1e-6 || q1->x > p1->x + 1e-6) {
	  fprintf (stderr, "ERROR: punto fuera de la arista\n");
	  fprintf (stderr, "       p1 = (%g, %g), p2 = (%g, %g)\n", 
		   p1->x, p1->y, p2->x, p2->y);
	  fprintf (stderr, "       q1 = (%g, %g), q2 = (%g, %g)\n", 
		   q1->x, q1->y, q2->x, q2->y);
	}
	if (q1->x < p2->x - 1e-6 || q2->x < p2->x - 1e-6) {
	  fprintf (stderr, "ERROR: punto fuera de la arista\n");
	  fprintf (stderr, "       p1 = (%g, %g), p2 = (%g, %g)\n", 
		   p1->x, p1->y, p2->x, p2->y);
	  fprintf (stderr, "       q1 = (%g, %g), q2 = (%g, %g)\n", 
		   q1->x, q1->y, q2->x, q2->y);
	}
      }
    }
    
    if (inter && (fabs (m * (q1->x - p1->x) + p1->y - q1->y) > 1e-6 ||
	fabs (m * (q2->x - p1->x) + p1->y - q2->y) > 1e-6)) {
      fprintf (stderr, "ERROR: punto no pertenece a linea\n");
      fprintf (stderr, "       p1 = (%g, %g), p2 = (%g, %g)\n", 
	       p1->x, p1->y, p2->x, p2->y);
      fprintf (stderr, "       q1 = (%g, %g), q2 = (%g, %g)\n", 
	       q1->x, q1->y, q2->x, q2->y);
    }
    
    if (inter && (q1->x > 1 || q1->x < 0 || q2->x > 1 || q2->x < 0 ||
		  q1->y > 1 || q1->y < 0 || q2->y > 1 || q2->y < 0)) {
      fprintf (stderr, "ERROR: punto fuera del cuadrado\n");
      fprintf (stderr, "       p1 = (%g, %g), p2 = (%g, %g)\n", 
	       p1->x, p1->y, p2->x, p2->y);
      fprintf (stderr, "       q1 = (%g, %g), q2 = (%g, %g)\n", 
	       q1->x, q1->y, q2->x, q2->y);
    }
    if (!inter) {
      fprintf (stderr, "NO\n");
    }
  }

  return 0;
}
