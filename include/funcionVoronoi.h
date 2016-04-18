#ifndef VIRPP_FUNCIONVORONOI
#define VIRPP_FUNCIONVORONOI

#include "funcion.h"
#include "mallaVoronoi.h"

class FuncionVoronoi: public virtual Funcion {
 protected:
  MallaVoronoi *malla;
 public:
  FuncionVoronoi (): Funcion (){}
  FuncionVoronoi (MallaVoronoi *m) {this->malla = m;}

  virtual void setMalla (MallaVoronoi *m) {
    if (malla != NULL) {
      eliminarMallaVoronoi (malla);
    }
    malla = m;
  }
  virtual MallaVoronoi *getMalla() {return malla;}

  virtual ~FuncionVoronoi () {
    if (malla != NULL) {
      //eliminarMallaVoronoi (malla);
    }
  }
};

#endif
