#ifndef VIRPP_INICIALIZADORPARAMETROS
#define VIRPP_INICIALIZADORPARAMETROS

#include "inicializador.h"
#include <iostream.h>
#include <cstdlib>
#include <fstream>
#include <string.h>

class InicializadorParametros: public Inicializador{
  double *pini;
 public:
  InicializadorParametros (double *p, int n);
  void inicializar (double pars[], int n, double init = 0.0);
  ~InicializadorParametros ();
};

#endif
