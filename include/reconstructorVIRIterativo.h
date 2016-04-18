#ifndef VIRPP_RECONSTRUCTORVIRITERATIVO
#define VIRPP_RECONSTRUCTORVIRITERATIVO

#include "reconstructorVIR.h"
#include "funcionVoronoiI.h"
#include "funcionVoronoiX.h"
#include "inicializadorParametros.h"

class ReconstructorVIRIterativo : public ReconstructorVIR {
 protected:
  /*char *nombreImagen, **nombresVis;
  int nVis, MEMIter, nPolsIni, nPolsFin, dn;
  double expanded_lx, expanded_ly, beamSize;*/
  int nIterMax;

  double realizarReconstruccion (FuncionVoronoi *funcCBI, double *p);

 public:
  ReconstructorVIRIterativo (char *nombreImagen, char **nombresVis, int nVis, 
			     int nPolsIni, int nPolsFin, int dn, 
			     double expanded_lx, double expanded_ly, int MEMIter = 1, 
			     double beamSize = -1, int nIterMax = 10);
  double run ();
};

#endif
