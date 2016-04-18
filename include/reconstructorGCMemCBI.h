#ifndef VIRPP_RECONSTRUCTORGCMEMCBI
#define VIRPP_RECONSTRUCTORGCMEMCBI

#include "reconstructor.h"
#include "reconstructorGC.h"
#include "funcionMemCBI.h"
#include "frprmn.h"

class ReconstructorGCMemCBI : public Reconstructor {
 protected:
  FuncionMemCBI *f;
  double ftol;
  int Nmax, nIterMax;
 public:
  ReconstructorGCMemCBI (FuncionMemCBI *f, Inicializador *ini, 
			 double ftol, int Nmax, double xini, int nIterMax = 1000000);
  ReconstructorGCMemCBI (char* nombe_archivo);
  double minimizador ();
  double run ();
  int getNmax() {return Nmax;};
  int getNIterMax() {return nIterMax;};
  FuncionMemCBI *getFunc() {return f;}
};

#endif
