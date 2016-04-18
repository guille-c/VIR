#ifndef VIRPP_INICIALIZADORARCHIVO
#define VIRPP_INICIALIZADORARCHIVO

#include "inicializador.h"
#include <iostream.h>
#include <cstdlib>
#include <fstream>
#include <string.h>

class InicializadorArchivo: public Inicializador{
  char nombre_archivo [255];
 public:
  InicializadorArchivo (char * nombre_archivo);
  void inicializar (double pars[], int n, double init = 0.0);
};

#endif
