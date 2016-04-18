#ifndef DISTRIBUCIONHAMILTONIANA
#define DISTRIBUCIONHAMILTONIANA

#include "distribucionPbb.h"
#include "distribucionGenerador.h"
#include "funcion.h"
#include <math.h>

class DistribucionHamiltoniana: public DistribucionGenerador {
 protected:
  DistribucionPbb *dMomentum;
  int Tau;
  double epsilon;
 public:
  DistribucionHamiltoniana (Funcion *P, DistribucionPbb *dMomentum,
			    int Tau, double epsilon);
  virtual double f (double pars[]) {return exp (-P->f(pars));}
  void generar (double *x_old);
  
  virtual ~DistribucionHamiltoniana () {delete dMomentum;}
};

#endif
