#ifndef DISTRIBUCIONMETROPOLIS
#define DISTRIBUCIONMETROPOLIS

#include "distribucionNormal.h"
#include "distribucionGenerador.h"
#include "funcion.h"
#include <math.h>

class DistribucionMetropolis: public DistribucionGenerador {
 protected:
  DistribucionNormal *Q;
 public:
  DistribucionMetropolis (Funcion *P, DistribucionNormal *Q);
  void generar (double *x_old);
  
  virtual ~DistribucionMetropolis () {delete Q;}
};

#endif
