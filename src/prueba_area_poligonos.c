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
//#include <gsl/gsl_multimin.h>

int main (int argc, char **argv) {
  int n = 1e3, i, j;
  float ftol = 1e-10;
  MallaVoronoi *malla;
  double sumaAreas = 0, A1, A2, x1, x2, x3, y1, y2, y3;
  NodoLista *nodo;
  PoligonoVoronoi *pol;
  double **im;
  double pols[1000][3];

  malla = newMallaVoronoi ();

    
  for (i = 0; i < 1e4; i++) {
    printf ("i: %d\n", i);
    eliminarMallaVoronoi (malla);
    malla = newMallaVoronoi ();
    for (j = 0; j < n; j++) {
      pols[j][0] = (double) random () / RAND_MAX;
      pols[j][1] = (double) random () / RAND_MAX;
      pols[j][2] = j;

      insertarSitio (malla, pols[j][0], pols[j][1], pols[j][2]);
      
      if (j != 0) {
	sumaAreas = 0;
	for (nodo = malla->poligonos->primero; nodo != NULL; nodo = nodo->sgte) {
	  pol = (PoligonoVoronoi *) nodo->info;
	  
	  sumaAreas += areaCuadrado (pol);
	}
	if (fabs(sumaAreas - 1) > 1e-10) {
	  fprintf (stderr, "ERROR suma == %g\n", sumaAreas);
	  fprintf (stderr, " diferencia = %g\n", sumaAreas - 1);
	  
	  for (nodo = malla->poligonos->primero; nodo != NULL; nodo = nodo->sgte) {
	    pol = (PoligonoVoronoi *) nodo->info;
	    
	    fprintf(stderr, "pol (%g, %g, %g) = %g\n", pol->x, pol->y, pol->valor, areaCuadrado (pol));
	  }
	  im = toImage(malla, 30, 30, NULL);
	  
	  for (j = 29; j >= 0; j--) {
	    for (i = 0; i < 30; i++)
	      {
		printf ("%g ", im[i][j]);
	      }
	    printf ("\n");
	  }
	  
	  exit(1);
	}
      }
    }
  }
  
  /* 
  for (nodo = malla->poligonos->primero; nodo != NULL; nodo = nodo->sgte) {
    pol = (PoligonoVoronoi *) nodo->info;
    
    sumaAreas += area (pol);
  }

  printf ("Suma areas = %g\n", sumaAreas);

  sumaAreas = 0;
  insertarSitio (malla, 0.156679, 0.400944, 12);
  insertarSitio (malla, 0.137232, 0.804177, 11);
  insertarSitio (malla, 0.0163006, 0.242887, 10);
  insertarSitio (malla, 0.141603, 0.606969, 9);
  insertarSitio (malla, 0.635712, 0.717297, 8);
  insertarSitio (malla, 0.95223, 0.916195, 7);
  insertarSitio (malla, 0.364784, 0.513401, 6);
  insertarSitio (malla, 0.477397, 0.628871, 5);
  insertarSitio (malla, 0.277775, 0.55397, 4);
  insertarSitio (malla, 0.335223, 0.76823, 3);
  insertarSitio (malla, 0.911647, 0.197551, 2);
  insertarSitio (malla, 0.783099, 0.79844, 1);
  insertarSitio (malla, 0.840188, 0.394383, 0);
  /*
  insertarSitio (malla, 0.5, 0.0, 1);
  insertarSitio (malla, 0.0, 1.0, 2);
  insertarSitio (malla, 1.0, 1.0, 3);
  
  insertarSitio (malla, 0.5, 0.5, 4);
    
  insertarSitio (malla, 0.2, 0.9, 5);
  insertarSitio (malla, 0.3, 0.4, 6);
  insertarSitio (malla, 0.2, 0.8, 7);
  insertarSitio (malla, 0.5, 0.4, 8);
  insertarSitio (malla, 1.0, 0.3, 9);
  insertarSitio (malla, 0.8, 0.2, 0);
  insertarSitio (malla, 0.6, 0.0, 1);
  insertarSitio (malla, 0.01, 1.0, 2);
  */
  for (nodo = malla->poligonos->primero; nodo != NULL; nodo = nodo->sgte) {
    double a;
    pol = (PoligonoVoronoi *) nodo->info;
    printf("\nsumando pol (%g, %g, %g) \n", pol->x, pol->y, pol->valor);
    a = areaCuadrado(pol);
    printf ("area: %g\n", a);
    sumaAreas += a;
  }

  printf ("Suma areas cuadrado = %g\n", sumaAreas);
  
  x1 = 0.5;
  y1 = -10;
  x2 = 10;
  y2 = 6.76419;
  x3 = -10;
  y3 = 7.46503;

  A1 = 0.5 * (x2 - x3 - (x3 - x1) / (y3 - y1) * (y2 - y3)) * (y3 - y2);
  printf ("Area 1: %g\n", A1);
  A2 = 0.5 * (x2 - x3 - (x3 - x1) / (y3 - y1) * (y2 - y3)) * (y2 - y1);
  printf ("Area 2: %g\n", A2);
  printf ("Area Real: %g\n", A1 + A2);
  
  im = toImage(malla, 20, 20, NULL);
  
  for (j = 19; j >= 0; j--) {
    for (i = 0; i < 20; i++)
      {
	printf ("%g ", im[i][j]);
      }
    printf ("\n");
  }
  
  printf("BORDE_INF = %d\n", BORDE_INF);
  printf("BORDE_DER = %d\n", BORDE_DER);
  printf("BORDE_SUP = %d\n", BORDE_SUP);
  printf("BORDE_IZQ = %d\n", BORDE_IZQ);
  /*
  Reconstructor *r;

  printf ("x = %g, %d\n", 1.23231, round (1.23231));
 
  r = newReconstructor (argv[1], argv + 2, argc - 2, n);
  //reinicializar (r, "MEM_NR_18.dat");
  run (r, 0.0, 0, ftol);
  */
  return 0;
}
