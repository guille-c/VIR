#ifndef INCREMENTADORMALLAARISTAS
#define INCREMENTADORMALLAARISTAS

#include "incrementador.h"
#include "mallaVoronoi.h"
#include "arista.h"

class IncrementadorMallaAristas: public Incrementador {
 protected:
  double valor_ini;
  
  void calcularNuevaPosicion(MallaVoronoi *m, double *x, double *y, double *valor);
 public:
  IncrementadorMallaAristas (double init_value);
  virtual double incrementar (double **p, Funcion *f);
};

#endif
