#ifndef VIRPP_FUNCIONCHI2VORONOIIMAGENN
#define VIRPP_FUNCIONCHI2VORONOIIMAGENN

#include <fstream>
#include <cstdlib>
#include <math.h>
#include "funcionVoronoi.h"
#include "funcionChi2VoronoiImagen.h"
#include "funcL.h"
#include "image_routines.h"
#include "mallaVoronoi.h"
#include "rasterizador.h"

class FuncionChi2VoronoiImagenN : public FuncionVoronoi {
 protected:
  FuncionChi2VoronoiImagen *func;
 public:
  FuncionChi2VoronoiImagenN (char * nombreImagen, int n) {
    func = new FuncionChi2VoronoiImagen (nombreImagen, n);
    this->n_pars = n * 3;
  }
  //FuncionChi2VoronoiImagen (char *nombre_archivo);
  double f (double pars[]) { return 1 / func->f(pars);}
  void df (double pars[], double grad[]) {func->df(pars, grad);}
  virtual void guardarInfo (string info) {func->guardarInfo (info);}

  virtual ~FuncionChi2VoronoiImagenN () {delete func;}

  //FuncionBayesVoronoiCBI& operator= (const FuncionBayesVoronoiCBI&);
};

#endif
