#ifndef VIRPP_FUNCIONPBBVORONOICBI
#define VIRPP_FUNCIONPBBVORONOICBI

#include "funcionCBI.h"
#include "funcionVoronoi.h"
#include "funcionBayesVoronoiCBI.h"

class FuncionPbbVoronoiCBI: public FuncionCBI, public FuncionVoronoi {
 protected:
  FuncionBayesVoronoiCBI *fBV;
  double c; // Constante de "normalizacion"
 public:
  FuncionPbbVoronoiCBI (char * nombreImagen, char **nombresVis, 
			int nVis, int n,
			int entropia, double cuantaSize,
			double expanded_lx, double expanded_ly);
  FuncionPbbVoronoiCBI (char *nombre_archivo);

  double f (double *pars);
  void df (double *pars, double * grad);
  void guardarInfo (char *info);
  void guardarInfo (string info);

  ~FuncionPbbVoronoiCBI();
};

#endif
