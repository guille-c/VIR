#ifndef VIRPP_FUNCIONCHI2VORONOIIMAGEN
#define VIRPP_FUNCIONCHI2VORONOIIMAGEN

#include <fstream>
#include <cstdlib>
#include <math.h>
#include "funcionVoronoi.h"
#include "funcL.h"
#include "image_routines.h"
#include "mallaVoronoi.h"
#include "rasterizador.h"

class FuncionChi2VoronoiImagen : public FuncionVoronoi {
 protected:
  struct image *im;
  char *nombreImagen;
 public:
  FuncionChi2VoronoiImagen (char * nombreImagen, int n_pols);
  //FuncionChi2VoronoiImagen (char *nombre_archivo);
  double f (double pars[]);
  void df (double pars[], double grad[]);
  virtual void guardarInfo (string info);

  virtual ~FuncionChi2VoronoiImagen ();

  //FuncionBayesVoronoiCBI& operator= (const FuncionBayesVoronoiCBI&);
};

#endif
