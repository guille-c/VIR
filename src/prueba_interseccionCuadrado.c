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
  int n = 3, i, j, nx = 128, ny = 128;
  float ftol = 1e-10;
  double **im;
  int **mask;
  MallaVoronoi *m;
  Reconstructor *r;

  //r = newReconstructor (argv[1], argv + 2, argc - 2, n);
  //reinicializar (r, "MEM_NR_065.dat");
  //run (r, 0.0, 0, ftol);

  PuntoVoronoi *p1 = newPuntoVoronoi (0.0, 1.3);
  PuntoVoronoi *p2 = newPuntoVoronoi (1.4, 0.5);
  PuntoVoronoi *q1 = newPuntoVoronoi (0.0, 1.1);
  PuntoVoronoi *q2 = newPuntoVoronoi (1.0, -0.1);
  double m1, b1;
  FILE *archivo = fopen("puntos.dat", "w");

  interseccionCuadrado (p1, p2, q1, q2, 4);
  
  printf ("p1 = (%g, %g)\n", p1->x, p1->y);
  printf ("p2 = (%g, %g)\n", p2->x, p2->y);
  m1 = (p1->y - p2->y) / (p1->x - p2->x);
  b1 = p1->y - m1 * p1->x;
  printf ("y = %gx + %g\n", m1, b1);
  printf ("q1 = (%g, %g)\n", q1->x, q1->y);
  printf ("q2 = (%g, %g)\n", q2->x, q2->y);

  fprintf (archivo, "%g %g\n", q1->x, q1->y);
  fprintf (archivo, "%g %g\n", q2->x, q2->y);

  return 0;
}
