#ifndef VIRPP_FUNCIONCERO
#define VIRPP_FUNCIONCERO

#include <iostream.h>
#include <fstream.h>
#include <string>
#include "funcion.h"
using namespace std;

class FuncionCero: public Funcion {
 protected:
 public:
  FuncionCero () {n_pars = 0;}
  FuncionCero (int n) {this->n_pars = n;}
  virtual double f(double pars[]) {return 0;}
  virtual void df(double pars[], double grad[]) {
    for (int i = 0; i < n_pars; i++) grad[i] = 0;
  }

  virtual ~FuncionCero() {}
};

#endif
