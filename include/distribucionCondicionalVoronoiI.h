#ifndef DISTRIBUCIONCONDICIONALVORONOII
#define DISTRIBUCIONCONDICIONALVORONOII

#include <math.h>
#include <string>
#include "distribucionCondicional.h"
#include "funcionBayesVoronoiCBI.h"
#include "inicializador.h"

class DistribucionCondicionalVoronoiI: public DistribucionCondicional {
 protected:
  DistribucionCondicional *dc;
  double *pos;
 public:
  DistribucionCondicionalVoronoiI (DistribucionCondicional *dc, Inicializador *ini);
  virtual double generar (double *, int i);
  virtual double evaluar (double *);
  virtual double evaluar (double *, double x, int i);

  virtual void guardarInfo (string s) {dc->guardarInfo(s);}

};

#endif
