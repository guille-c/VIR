#ifndef VIRPP_FUNCIONEXP
#define VIRPP_FUNCIONEXP

#include <math.h>
#include "funcion.h"

class FuncionExp: public Funcion {
 public:
  FuncionExp (int n): Funcion (n) {}
  virtual double f(double pars[]);
  virtual void df(double pars[], double grad[]);
};

#endif
