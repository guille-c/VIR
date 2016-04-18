#ifndef VIRPP_FUNCIONMEMCBI
#define VIRPP_FUNCIONMEMCBI

#include <fstream>
#include <cstdlib>
#include <math.h>
#include "funcionCBI.h"
#include "image_routines.h"
#include "uvsubs.h"
#include "mockcbiRoutines.h"

class FuncionMemCBI : public FuncionCBI {
 protected:
  Funcion *entropia;
  double chi2, S;
  //double minpix;
  double lambda, pMin, noise_cut;

 public:
  FuncionMemCBI (char *nombreImagen, char **nombresVis,
		 int nVis, Funcion *entropia, double Imin = 0.0,
		 double lambda = 1.0, double noise_cut = 1e100, double beamSize = -1);
  //FuncionMemCBI (char *nombre_archivo);
  double f (double pars[]);
  void df (double pars[], double grad[]);
  virtual void guardarInfo (string info);
  //virtual void setQuanta (double quantaSize) {ruido = quantaSize;}
  //virtual double getQuanta () {return ruido;}
  double getChi2 () {return chi2;}
  double getS () {return S;}
  double getLambda () {return lambda;}
  void setLambda (double lambda);
  virtual void resizeImage (int nx, int ny, int cx, int cy);
  void scalePars (double *pars);

  FuncionMemCBI& operator= (const FuncionMemCBI&);
};

#endif
