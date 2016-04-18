#ifndef VIRPP_FUNCIONVORONOII
#define VIRPP_FUNCIONVORONOII

#include "funcion.h"
#include "funcionVoronoi.h"
#include "inicializador.h"

class FuncionVoronoiI: public Funcion {
 protected:
  double *pos;
  FuncionVoronoi *fV;

 public:
  FuncionVoronoiI (FuncionVoronoi *f, Inicializador *ini);
  FuncionVoronoiI (FuncionVoronoi *f, double *pos);

  virtual double f (double *pars);
  virtual void df (double *pars, double * grad);
  virtual void guardarInfo (string info);
  virtual void setNPars (int n_pars);

  virtual ~FuncionVoronoiI();
};

#endif
