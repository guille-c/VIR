#include <iostream.h>
#include <stdlib.h>
#include <fstream.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "image_routines.h"
#include "rasterizador.h"

void actualizarIntensidades (struct image *im, double *p, int nPols);

int main(int argc, char **argv) {
  std::ifstream archivoSitios (argv[1]);
  struct image *im = do_read(argv[2]);
  std::ofstream archivoMalla (argv[3]);
  double x, y;  
  int nPols, i, j, nx = im->size[0], ny = im->size[1];
  
  archivoSitios >> nPols;
  cout << "nPols = " << nPols << "\n";
  
  double *p = new double [3 * nPols];
  for (int cont = 0; cont < nPols; cont++) {
    archivoSitios >> i;
    archivoSitios >> j;
    p[3*cont    ] = i / (nx - 1.0);
    p[3*cont + 1] = j / (ny - 1.0);
    p[3*cont + 2] = 0;
  }
  actualizarIntensidades (im, p, nPols);

  archivoMalla << nPols << "\n";
  MallaVoronoi *m = newMallaVoronoi ();
  for (int cont = 0; cont < nPols; cont++) {
    insertarSitio (m, &p[3*cont], &p[3*cont+1], p[3*cont+2]);
    archivoMalla << p[3*cont    ] << "\t";
    archivoMalla << p[3*cont + 1] << "\t";
    archivoMalla << p[3*cont + 2] << "\n";
  }
  
  double **im2 = toImage (m, nx, ny, NULL);
  
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      im->pixels[i + nx * j] = im2[i][j];
    }
  }
  do_write_fits (im, argv[4]);
}

void actualizarIntensidades (struct image *im, double *p, int nPols) {
  MallaVoronoi *malla = newMallaVoronoi();
  int i, j, nx = im->size[0], ny = im->size[1], *conts;
  double valor, *valores, x, y;
  PoligonoVoronoi *pol;

  // Creamos la malla antigua
  for (i = 0; i < nPols; i++) {
    insertarSitio (malla, &(p[3 * i]), &(p[3 * i + 1]), p[3 * i + 2]);
  }

  valores = new double [nPols];
  conts = new int [nPols];
  for (i = 0; i < nPols; i++) {
    valores[i] = 0;
    conts[i] = 0;
  }
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      x = i / (nx - 1.0);
      y = j / (ny - 1.0);

      pol = encontrarPoligono (malla, x, y);
      valores [pol->id - 3]  += im->pixels[i + j*nx];
      conts [pol->id - 3] ++;
      }
    }
  
  eliminarMallaVoronoi(malla);  
  for (i = 0; i < nPols; i++) {
    if (conts[i] != 0) {
      p[3*i + 2] = valores[i] / conts[i];
    }
    else {
      p[3*i + 2] = 0;
    }
  }
  delete [] valores;
  delete [] conts;
}
