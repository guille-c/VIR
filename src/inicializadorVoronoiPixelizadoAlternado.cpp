#include "inicializadorVoronoiPixelizadoAlternado.h"

void InicializadorVoronoiPixelizadoAlternado::inicializar (double pars[], int n,
							   double init){
  int i, j, k = 0, n_pols, nx, ny, n_cajas, parx, pary;
  double dx, dy, x, y, delta = 1e-6;

  n_pols = n / 3;
  n_cajas = n_pols / 10;
  nx = (int) round (sqrt (n_cajas));
  ny = (int) round (sqrt (n_cajas));
  if (nx * ny * 3 * 10 != n) {
    cerr << "ERROR en InicializadorVoronoiPixelizado::inicializar, nx * ny * 30 != n\n";
    cerr << "      " << nx << "*" << ny << "*30 != " << n << "\n";
    exit (1);
  }
  
  // primero ponemos los poligonos grandes
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny * 2; j++) {
      x = (2*(j % 2) + 1.0 + i * 4) / (4.0 * nx);
      y = (1.0 + j * 2) / (4.0 * ny);
      pars[3 * k]     = x;
      pars[3 * k + 1] = y;
      pars[3 * k + 2] = init;
      k++;
    }
  }

  // ahora, los chicos.
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny * 2; j++) {
      x = (2*((j + 1) % 2) + 1.0 + i * 4) / (4.0 * nx);
      y = (1.0 + j * 2) / (4.0 * ny);

      pars[3 * k]     = x + delta;
      pars[3 * k + 1] = y + delta;
      pars[3 * k + 2] = init;
      k++;

      pars[3 * k]     = x - delta;
      pars[3 * k + 1] = y + delta;
      pars[3 * k + 2] = init;
      k++;

      pars[3 * k]     = x + delta;
      pars[3 * k + 1] = y - delta;
      pars[3 * k + 2] = init;
      k++;

      pars[3 * k]     = x - delta;
      pars[3 * k + 1] = y - delta;
      pars[3 * k + 2] = init;
      k++;
    }
  }
}
