#ifndef DISTRIBUCIONCONDICIONALMEMCBI
#define DISTRIBUCIONCONDICIONALMEMCBI

#include <math.h>
#include <string>
#include "distribucionCondicionalCBI.h"
#include "funcion.h"
#include "funcionMemCBI.h"
#include "inicializador.h"

class DistribucionCondicionalMemCBI: public DistribucionCondicionalCBI {
 protected:
  Funcion *entropia;
  double chi2, S;
  double minpix;
  double lambda;
  double L_old;
  int N;

  double calcularChi2();
 public:
  DistribucionCondicionalMemCBI (char *nombreImagen, char **nombresVis, int nVis, 
				 Funcion *entropia, double lambda, int N, Inicializador *ini = 0);
  virtual double generar (double *, int i);
  virtual double evaluar (double *);
  virtual double evaluar (double *, double x, int i);
  virtual double evaluar_old (double *, double x, int i);

  //virtual void guardarInfo (string s) {f->guardarInfo(s);}

  virtual ~DistribucionCondicionalMemCBI();
};

#endif
