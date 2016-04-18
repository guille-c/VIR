#ifndef DISTRIBUCIONCONDICIONALBAYESVORONOICBII
#define DISTRIBUCIONCONDICIONALBAYESVORONOICBII

#include <math.h>
#include <string>
#include "distribucionCondicionalBVCBI.h"
#include "funcionBayesVoronoiCBI.h"
#include "inicializador.h"

class DistribucionCondicionalBVCBI_I: public DistribucionCondicionalBVCBI {
 protected:
  double *pos;
 public:
  DistribucionCondicionalBVCBI_I (char *nombreArchivoInput, int N, Inicializador *ini);
  virtual double generar (double *, int i);
  virtual double evaluar (double *, double x, int i);

};

#endif
