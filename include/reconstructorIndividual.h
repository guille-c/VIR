#ifndef VIRPP_RECONSTRUCTORINDIVIDUAL
#define VIRPP_RECONSTRUCTORINDIVIDUAL

#include <math.h>
#include <iostream.h>
#include <fstream.h>
#include <stdarg.h>
#include <string>
#include "reconstructor.h"
#include "inicializador.h"
#include "funcion.h"

class ReconstructorIndividual: public Reconstructor {
 protected:
  double *max, *min;
  double *N;
  Funcion *f;
  //Lo de Reconstructor
  //int iter, n_pars;
  //double *p;
  //Inicializador *ini;

 public:
  ReconstructorIndividual (Funcion *f, double min, double max, double N, Inicializador *ini);
  ReconstructorIndividual (Funcion *f, double *min, double *max, double *N, Inicializador *ini);
  virtual double run ();

  //void setNPars (double n_pars) {this->n_pars = n_pars;}

  virtual ~ReconstructorIndividual();
  
};

#endif
