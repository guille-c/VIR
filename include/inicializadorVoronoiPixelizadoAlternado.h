#ifndef VIRPP_INICIALIZADORVORONOIPIXELIZADOALTERNADO
#define VIRPP_INICIALIZADORVORONOIPIXELIZADOALTERNADO

#include "inicializador.h"
#include <cstdlib>
#include <math.h>
#include <iostream.h>

class InicializadorVoronoiPixelizadoAlternado: public Inicializador{
 public:
  void inicializar (double pars[], int n, double init = 0.0);
};

#endif
