#ifndef INCREMENTADORMALLAPIXEL
#define INCREMENTADORMALLAPIXEL

#include "incrementador.h"
#include "funcionBayesVoronoiCBI.h"
#include "mallaVoronoi.h"

class IncrementadorMallaPixel: public Incrementador {
 protected:
  int nx, ny;
  double valor_ini;
 public:
  IncrementadorMallaPixel (int nx, int ny, double init_value);
  virtual double incrementar (double **p, Funcion *f);
};

#endif
