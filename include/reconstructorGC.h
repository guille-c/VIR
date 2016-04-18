#ifndef VIRPP_RECONSTRUCTORGC
#define VIRPP_RECONSTRUCTORGC

#include "reconstructor.h"
#include "funcion.h"
#include "frprmn.h"

class ReconstructorGC : public Reconstructor {
 protected:
  Funcion *f;
  double ftol;
 public:
  ReconstructorGC (Funcion *f, Inicializador *ini, double ftol);
  ReconstructorGC (char* nombe_archivo);
  double minimizador ();
  double run ();

  virtual void setNPars (int n_pars);
  virtual void guardarInfo (string info) {f->guardarInfo(info);}

  Funcion * getFunc() {return f;}
};

#endif
