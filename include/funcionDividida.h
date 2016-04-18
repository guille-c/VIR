#ifndef VIRPP_FUNCIONDIVIDIDA
#define VIRPP_FUNCIONDIVIDIDA

#include <iostream.h>
#include <fstream.h>
#include <string>
#include <math.h>
#include "funcion.h"

class FuncionDividida: public Funcion {
 protected:
  Funcion *funcion;
 public:
  FuncionDividida (Funcion *f);
  virtual double f(double pars[]);
  virtual void df(double pars[], double grad[]);
  virtual void guardarInfo (string info);
  virtual void setNPars (int n_pars);

  virtual ~FuncionDividida();
};

#endif
