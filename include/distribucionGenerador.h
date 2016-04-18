#ifndef DISTRIBUCIONGENERADOR
#define DISTRIBUCIONGENERADOR

#include "distribucionPbb.h"
#include "funcion.h"
#include <math.h>

class DistribucionGenerador: public DistribucionPbb {
 protected:
  Funcion *P;
 public:
  DistribucionGenerador (Funcion *P): DistribucionPbb(P->getNPars()) {
    this->P = P;
  }
  virtual double f (double pars[]) {return P->f(pars);}
  void setFuncion (Funcion *P) {this->P = P; this->n_pars = P->getNPars();}
  Funcion *getFuncion () {return P;}
  virtual void guardarInfo(string s) {P->guardarInfo(s);}

  virtual ~DistribucionGenerador () {}//{delete P; P = NULL;}
};

#endif
