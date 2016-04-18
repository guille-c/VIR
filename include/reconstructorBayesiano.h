#ifndef RECONSTRUCTORBAYESIANO
#define RECONSTRUCTORBAYESIANO

#include "reconstructor.h"
#include "distribucionGenerador.h"

class ReconstructorBayesiano: public Reconstructor {
 protected:
  DistribucionPbb *distribucion;
  int n_iter;
 public:
  ReconstructorBayesiano (DistribucionPbb *d, Inicializador *ini, int n_iter);
  double run ();

  ~ReconstructorBayesiano();
};

#endif
