#ifndef VIRPP_FUNCIONRECONSTRUCTOR
#define VIRPP_FUNCIONRECONSTRUCTOR

#include "funcion.h"
#include "reconstructor.h"

using namespace std;

class FuncionReconstructor: public Funcion {
 protected:
  Reconstructor *rec;
 public:
  FuncionReconstructor () : Funcion() {}
  FuncionReconstructor (Reconstructor *rec);
  virtual double f(double pars[]);
  virtual void df(double pars[], double grad[]);
  virtual void guardarInfo (string info);
  virtual void setNPars(int n_pars);

  virtual ~FuncionReconstructor();
};

#endif
