#ifndef VIRPP_FUNCIONBAYESVORONOICBII
#define VIRPP_FUNCIONBAYESVORONOICBII

#include <fstream>
#include <cstdlib>
#include <math.h>
#include "funcionCBI.h"
#include "funcionVoronoi.h"
#include "inicializador.h"
#include "funcL.h"

class FuncionBayesVoronoiCBII : public FuncionCBI, public FuncionVoronoi {
 protected:
  double *pos;
  funcL *fL;
 public:
  FuncionBayesVoronoiCBII (char * nombreImagen, char **nombresVis,
			   int nVis, int n,
			   double init_value, int init_gauss, 
			   int entropia, double cuantaSize,
			   double expanded_lx, double expanded_ly, double beamSize = -1);
  FuncionBayesVoronoiCBII (char *nombre_archivo, Inicializador *ini);
  double f (double pars[]);
  void df (double pars[], double grad[]);
  virtual void guardarInfo (string info);
  virtual void setNPars (int n_pars);

  MallaVoronoi *getMalla() {return fL->malla;}

  virtual ~FuncionBayesVoronoiCBII ();

  FuncionBayesVoronoiCBII& operator= (const FuncionBayesVoronoiCBII&);
};

#endif
