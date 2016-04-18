#ifndef VIRPP_INICIALIZADORVORONOIPIXELIZADO
#define VIRPP_INICIALIZADORVORONOIPIXELIZADO

#include "inicializador.h"
#include <cstdlib>
#include <math.h>
#include <iostream.h>

class InicializadorVoronoiPixelizado: public Inicializador{
 public:
  void inicializar (double pars[], int n, double init = 0.0);
};

#endif
