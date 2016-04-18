#include "inicializadorVoronoiUniforme.h"

void InicializadorVoronoiUniforme::inicializar (double pars[], int n,
						double init){
  int i, n_pols;

  n_pols = n / 3;
  if (pars == NULL) {
    pars = new double [n];
  }
  for (i = 0; i < n_pols; i++) {
    pars[3 * i]  = rand()/(RAND_MAX + 1.0);
    pars[3 * i + 1]  = rand()/(RAND_MAX + 1.0);
    pars[3 * i + 2]  = init;
  }
}
