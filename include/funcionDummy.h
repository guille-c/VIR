#ifndef VIRPP_FUNCIONDUMMY
#define VIRPP_FUNCIONDUMMY

#include <iostream.h>
#include <fstream.h>
#include "funcion.h"

class FuncionDummy: public Funcion {
 public:
  double f(double pars[]) {
    double ret;
    for (int i = 0; i < n_pars; i++) {
      ret += pars[i] * pars[i];
    }
    return ret;
  }
  void df(double pars[], double grad[]) {
    double ret;
    for (int i = 0; i < n_pars; i++) {
      grad[i] = 2 * pars[i];
    }
  }
};

#endif
