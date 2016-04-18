#ifndef VIRPP_INICIALIZADORCONSTANTE
#define VIRPP_INICIALIZADORCONSTANTE

#include "inicializador.h"

class InicializadorConstante: public Inicializador {
 public:
  virtual void inicializar (double pars[], int n, double init = 0.0) {
    for (int i = 0; i < n; i++) {
      pars[i] = init;
    }
  }
};

#endif
