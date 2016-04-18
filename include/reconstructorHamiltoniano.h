#ifndef VIRPP_RECONSTRUCTORHAMILTONIANO
#define VIRPP_RECONSTRUCTORHAMILTONIANO

#include "reconstructor.h"
#include "distribucionGenerador.h"

class ReconstructorHamiltoniano: public Reconstructor {
 protected:
  DistribucionGenerador *distribucionMomentum;
  int n_iter, Tau;
  double epsilon;
 public:
  ReconstructorHamiltoniano (Funcion *f, DistribucionGenerador *d, Inicializador *ini,
			     int n_iter, int Tau, double epsilon);
  double run ();

  ~ReconstructorHamiltoniano();
};

#endif
