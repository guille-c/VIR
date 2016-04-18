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

  mask = (int **) malloc (nx * sizeof (int *));
  for (i = 0; i < 128; i++) {
    mask [i] = (int *) malloc (ny * sizeof (int));
  }
  
  for (i = 0; i < 1e3; i++) {
    printf ("i = %d\n", i);
    m = newMallaVoronoi();
    insertarSitio(m, (double) rand() / RAND_MAX,  
		  (double) rand() / RAND_MAX, j);
    for (j = 0; j < 1e2; j++) {
      insertarSitio(m, (double) rand() / RAND_MAX,  
		    (double) rand() / RAND_MAX, j);
      im = toImage(m, nx, ny, mask);
      if (!comprobarImagen (m, mask, nx, ny)){
	guardarFits (im, mask, "ERRORcomprobarImagen", argv[1]);
	im = toImageSinRaster(m, nx, ny, mask);
	guardarFits (im, mask, "ERRORcomprobarImagenSinRaster", argv[1]);
	imprimirMallaArchivo (m, "ERRORcomprobarImagen.dat");
	exit (1);
      }
      for (j = 0; j < nx; j++) {
	free (im[j]);
      }
      free (im);
    }
    if (m != NULL) {
      eliminarMallaVoronoi (m);
      m = NULL;
    }
  }

  return 0;
}
