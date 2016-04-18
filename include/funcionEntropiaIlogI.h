#ifndef VIRPP_FUNCIONENTROPIAILOGI
#define VIRPP_FUNCIONENTROPIAILOGI

#include <iostream.h>
#include <fstream.h>
#include <string>
#include <math.h>
#include "funcion.h"
using namespace std;

class FuncionEntropiaIlogI: public Funcion {
  double M, lambda;
 public:
  FuncionEntropiaIlogI (int n, double M, double l) {this->n_pars = n; this->M = M; lambda = l;}
  virtual double f(double pars[]);
  virtual void df(double pars[], double grad[]);
};

#endif
