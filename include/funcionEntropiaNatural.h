#ifndef VIRPP_FUNCIONENTROPIANATURAL
#define VIRPP_FUNCIONENTROPIANATURAL

#include <iostream.h>
#include <fstream.h>
#include <string>
#include <math.h>
#include "funcion.h"
using namespace std;

class FuncionEntropiaNatural: public Funcion {
 public:
  FuncionEntropiaNatural (int n, double lambda = 1) {this->n_pars = n;}
  virtual double f(double pars[]);
  virtual void df(double pars[], double grad[]);
};

#endif
