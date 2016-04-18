#ifndef VIRPP_FUNCIONNEGATIVA
#define VIRPP_FUNCIONNEGATIVA

#include <iostream.h>
#include <fstream.h>
#include <string>
#include "funcion.h"

class FuncionNegativa: public Funcion {
 protected:
  Funcion *funcion;
 public:
  FuncionNegativa (Funcion *f);
  virtual double f(double pars[]);
  virtual void df(double pars[], double grad[]);
  virtual void guardarInfo (string info);
  virtual void setNPars (int n_pars);

  virtual ~FuncionNegativa();
};

#endif
