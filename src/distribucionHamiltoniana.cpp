#include "distribucionHamiltoniana.h"

DistribucionHamiltoniana::DistribucionHamiltoniana (Funcion *E, DistribucionPbb *dMomentum,
						    int Tau, double epsilon)
  : DistribucionGenerador (E) {
  if (P->getNPars() != dMomentum->getNPars()) {
    cerr << "ERROR: P y Q tienen distinto numero de parametros\n";
    cerr << "       " << P->getNPars() << " vs " << dMomentum->getNPars() << "\n";
    exit (1);
  }
  this->dMomentum = dMomentum;
  this->Tau = Tau;
  this->epsilon = epsilon;
}

void DistribucionHamiltoniana::generar (double *x_old) {
  double *x_new = new double [n_pars];
  double *p = new double [n_pars];
  double *g = new double [n_pars];
  double *g_new = new double [n_pars];
  double E, H;
 
  E = P->f(x_old);
  P->df (x_old, g);

  dMomentum->generar(p);

  for (int i = 0; i < n_pars; i++) {
    x_new[i] = x_old[i];
    g_new[i] = g[i];
    H += p[i] * p[i] / 2;
  }
  H += E;

  for (int tau = 0; tau < Tau; tau++) {
    cout << "--------------tau = " << tau << " --------\n";
    for (int i = 0; i < n_pars; i++) {
      p[i] -= epsilon * g_new[i] / 2;
      x_new[i] += epsilon * p[i];
    }

    P->df (x_new, g_new);

    for (int i = 0; i < n_pars; i++) {
      p[i] -= epsilon * g_new[i] / 2;
    }    
  }

  double E_new = P->f (x_new);
  double H_new = E_new;
  for (int i = 0; i < n_pars; i++) {
    H_new += p[i] * p[i] / 2;
  }
  double dH = H_new - H;
  
  if (dH < 0 || (double) rand () / RAND_MAX < exp (-dH)) {
    for (int i = 0; i < n_pars; i++) {
      x_old[i] = x_new[i];
    }
  }
  delete [] x_new;
  delete [] g;
  delete [] g_new;
}
