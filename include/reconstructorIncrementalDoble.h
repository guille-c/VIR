#ifndef RECONSTRUCTORINCREMENTALDOBLE
#define RECONSTRUCTORINCREMENTALDOBLE

#include "reconstructor.h"
#include "reconstructorIncremental.h"
#include "incrementador.h"

class ReconstructorIncrementalDoble: public ReconstructorIncremental {
 protected:
  //Funcion *f;
  //Incrementador *inc;
  //int n_iter;
  Reconstructor *r;
 public:
  ReconstructorIncrementalDoble (Funcion *f, Inicializador *ini, 
				 Incrementador *inc, int n_iter, Reconstructor *r);
  virtual double run ();

  virtual ~ReconstructorIncrementalDoble();
};

#endif
