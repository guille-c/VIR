#ifndef INCREMENTADORMALLAIMAGEN
#define INCREMENTADORMALLAIMAGEN

#include "incrementador.h"
#include "mallaVoronoi.h"
#include "image_routines.h"

class IncrementadorMallaImagen: public Incrementador {
 protected:
  struct image *im;
  double valor_ini;
  void encontrarMayorError (struct image *im, MallaVoronoi *m, 
			    int *imax, int *jmax, double *valor);
  void encontrarMayorError2 (struct image *im, MallaVoronoi *m, 
			    int *imax, int *jmax, double *valor);
  void ActualizarIntensidades (double *p, int nPols);
 public:
  IncrementadorMallaImagen (char * nombreImagen, double init_value);
  virtual double incrementar (double **p, Funcion *f);

  virtual ~IncrementadorMallaImagen();
};

#endif
