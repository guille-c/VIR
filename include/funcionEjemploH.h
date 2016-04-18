#ifndef VIRPP_FUNCIONEJEMPLOH
#define VIRPP_FUNCIONEJEMPLOH

#include <math.h>
#include "funcion.h"

class FuncionEjemploH: public Funcion {
 public:
  FuncionEjemploH (int n): Funcion (n) {}
  virtual double f(double pars[]);
  virtual void df(double pars[], double grad[]);
};

#endif
