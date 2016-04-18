#include "distribucionMetropolis.h"

DistribucionMetropolis::DistribucionMetropolis (Funcion *P, DistribucionNormal *Q)
  : DistribucionGenerador (P) {
  if (P->getNPars() != Q->getNPars()) {
    cerr << "ERROR: P y Q tienen distinto numero de parametros\n";
    cerr << "       " << P->getNPars() << " vs " << Q->getNPars() << "\n";
    exit (1);
  }
  this->Q = Q;
}

void DistribucionMetropolis::generar (double *x_old) {
  double *x_new = new double [n_pars];
  double a;
 
  Q->setMedia(x_old);
  Q->generar(x_new);
  a = P->f(x_new)/ P->f(x_old);
  if (a >= 1 || (double) rand () / RAND_MAX <= a) { //aceptado
    for (int i = 0; i < n_pars; i++) {
      x_old[i] = x_new[i];
    }
  }
  delete [] x_new;
}
