#ifndef RECONSTRUCTORINCREMENTAL
#define RECONSTRUCTORINCREMENTAL

#include "reconstructor.h"
#include "incrementador.h"

class ReconstructorIncremental: public Reconstructor {
 protected:
  Funcion *f;
  Incrementador *inc;
  int n_iter;
 public:
  ReconstructorIncremental (Funcion *f, Inicializador *ini, 
			    Incrementador *inc, int n_iter);
  virtual double run ();

  virtual ~ReconstructorIncremental();
};

#endif
