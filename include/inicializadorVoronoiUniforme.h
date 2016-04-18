#ifndef VIRPP_INICIALIZADORVORONOIUNIFORME
#define VIRPP_INICIALIZADORVORONOIUNIFORME

#include "inicializador.h"
#include <cstdlib>

class InicializadorVoronoiUniforme: public Inicializador{
 public:
  void inicializar (double pars[], int n, double init = 0.0);
};

#endif
