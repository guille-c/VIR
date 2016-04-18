#ifndef DISTRIBUCIONGIBBS
#define DISTRIBUCIONGIBBS

#include "distribucionCondicional.h"
#include "distribucionPbb.h"
#include "funcion.h"
//#include "reconstructor.h"  // solo para i2s
#include <math.h>
#include <string>

class DistribucionGibbs: public DistribucionPbb {
 protected:
  DistribucionCondicional *Pcond;
 public:
  DistribucionGibbs (DistribucionCondicional *Pcond);
  virtual double f (double *);
  virtual void generar (double *x);
  virtual void guardarInfo (string s) {Pcond->guardarInfo(s);}
};

#endif
