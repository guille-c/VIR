#include "reconstructorHamiltoniano.h"

ReconstructorHamiltoniano::ReconstructorHamiltoniano (Funcion *f, DistribucionGenerador *d, 
						      Inicializador *ini, int n_iter, 
						      int Tau, double epsilon)
  : Reconstructor (f, ini) {
  this->distribucionMomentum = d;
  this->n_iter = n_iter;
  this->Tau = Tau;
  this->epsilon = epsilon;
}

double ReconstructorHamiltoniano::run () {
  double f_mejor = -1e300, func;
  double *p_aux = new double [f->getNPars()];
  double *g = new double [f->getNPars()];

  for (int j = 0; j < f->getNPars(); j++) {
    p_aux[j] = p[j];
    cout << "p_aux[" << j << "] = " << p_aux[j] << "\n";
  }

  for (iter = 0; iter < n_iter; iter++) {
    if ((func = distribucion->f(p_aux)) > f_mejor) {
      cout << "Mejor en iteracion " << iter << ", " << p_aux[0] << "\n";
      cout << "valor: " << func << "\n";
      f->guardarInfo (i2s (iter, 9));
      f_mejor = func;
      imprimirLog ("MEJOR: iteracion %d: %g en\n", iter, f_mejor);
      for (int j = 0; j < distribucion->getNPars(); j++) {
	imprimirLog ("  p[%d] = %g\n", j, p_aux[j]);
	p[j] = p_aux[j];
      }
    }
    //cout << func << " as funcion\n";
    distribucion->generar (p_aux);

    //imprimimos en reconstructor.log
    if (iter % 100 == 0) {
      cout << "iter: " << iter << "\n";
      imprimirLog ("Iteracion %d, calculando funcion con parametros:\n", iter);
      for (int j = 0; j < distribucion->getNPars(); j++) {
	imprimirLog ("  p[%d] = %g\n", j, p_aux[j]);
      }
      imprimirLog (" valor funcion: %g, funcion mejor: %g\n", func, f_mejor);
      }
  }
  delete [] p_aux;
  return 0.0;
}

ReconstructorHamiltoniano::~ReconstructorHamiltoniano() {
  delete distribucionMomentum;
}
