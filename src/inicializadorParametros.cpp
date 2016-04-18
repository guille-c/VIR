#include "inicializadorParametros.h"

InicializadorParametros::InicializadorParametros (double *p, int n) {
  pini = new double [n];
  for (int i = 0; i < n; i++) {
    pini [i] = p [i];
  }
}

void InicializadorParametros::inicializar (double pars[], int n, double init) {
  for (int i = 0; i < n; i++) {
    pars [i] = pini[i];
  }
  
}

InicializadorParametros::~InicializadorParametros () {
  delete [] pini;
}
