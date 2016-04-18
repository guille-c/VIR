#ifndef VIRPP_FUNCIONVORONOIX
#define VIRPP_FUNCIONVORONOIX

#include "funcion.h"
#include "funcionVoronoi.h"
#include "inicializador.h"

class FuncionVoronoiX: public Funcion {
 protected:
  double *I;
  FuncionVoronoi *fV;

 public:
  FuncionVoronoiX (FuncionVoronoi *f, Inicializador *ini);
  FuncionVoronoiX (FuncionVoronoi *f, double *I);

  virtual double f (double *pars);
  virtual void df (double *pars, double * grad);
  virtual void guardarInfo (string info);
  virtual void setNPars (int n_pars);

  virtual ~FuncionVoronoiX();
};

#endif
