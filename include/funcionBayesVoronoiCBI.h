#ifndef VIRPP_FUNCIONBAYESVORONOICBI
#define VIRPP_FUNCIONBAYESVORONOICBI

#include <fstream>
#include <cstdlib>
#include <math.h>
#include "funcionCBI.h"
#include "funcionVoronoi.h"
#include "funcL.h"

class FuncionBayesVoronoiCBI : public FuncionCBI, public FuncionVoronoi {
 protected:
  funcL *fL;
 public:
  FuncionBayesVoronoiCBI (char * nombreImagen, char **nombresVis, 
			  int nVis, int n,
			  int entropia, double cuantaSize,
			  double lx, double ly, double beamSize = -1);
  FuncionBayesVoronoiCBI (char *nombre_archivo);
  double f (double pars[]);
  void df (double pars[], double grad[]);
  virtual void guardarAtenuaciones ();
  virtual void guardarInfo (string info);
  virtual void setNPars (int n_pars);
  virtual void setRuido (double quantaSize) {fL->difmapNoise = quantaSize; ruido = quantaSize;}
  virtual double getRuido () {return fL->difmapNoise;}
  funcL *getFL() {return fL;}

  MallaVoronoi *getMalla() {return fL->malla;}

  virtual ~FuncionBayesVoronoiCBI ();

  FuncionBayesVoronoiCBI& operator= (const FuncionBayesVoronoiCBI&);
};

#endif
