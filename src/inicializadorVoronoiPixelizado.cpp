#include "inicializadorVoronoiPixelizado.h"

void InicializadorVoronoiPixelizado::inicializar (double pars[], int n,
						  double init){
  int i, j, n_pols, nx, ny;

  n_pols = n / 3;
  nx = (int) sqrt (n_pols);
  ny = (int) sqrt (n_pols);
  if (nx * ny * 3 != n) {
    cerr << "ERROR en InicializadorVoronoiPixelizado::inicializar, nx * ny * 3 != n\n";
    exit (1);
  }
  if (pars == NULL) {
    pars = new double [3 * n_pols];
  }
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      pars[3 * (i + j * nx)]     = i / (nx - 1.0);
      pars[3 * (i + j * nx) + 1] = j / (ny - 1.0);;
      pars[3 * (i + j * nx) + 2] = init;
    }
  }
}
