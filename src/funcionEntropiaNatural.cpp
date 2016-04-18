#include "funcionEntropiaNatural.h"
/*
double FuncionEntropiaNatural::f(double pars[]) {
  int i;
  double S = 0, N = 0, Ni, pMin = 1e300;
  
  for (i = 0; i < n_pars; i++) {
    if (pars[i] < pMin) {
      pMin = pars[i];
    }
  }

  for (i = 0; i < n_pars; i++) {
    Ni = pars[i] - pMin;
    //N += pars[i];
    N += Ni;
    //S -= lgamma(pars[i] + 1);
    S -= lgamma(Ni + 1);
  }
  S += lgamma (N + 1) - N * log ((double) n_pars);

  if (isinf(S)) {
    S = 0;
  }
  
  return S;
}


void FuncionEntropiaNatural::df(double pars[], double grad[]){
  int i, j, imin;
  double dS_cte = 0, N = 0, Ni, pMin = 1e300;
  
  for (i = 0; i < n_pars; i++) {
    if (pars[i] < pMin) {
      pMin = pars[i];
      imin = i;
    }
  }

  for (i = 0; i < n_pars; i++) {
    Ni = pars [i] - pMin;
    N += Ni;
    //N += pars[i];
  }

  dS_cte = -log ((double) n_pars);
  for (i = 1; i <= N; i++) {
    dS_cte += 1.0 / i;
  }
  
  for (i = 0; i < n_pars; i++) {
    grad[i] = dS_cte;
    Ni = pars [i] - pMin;
    //for (j = 1; j <= pars[i]; j++) {
    for (j = 1; j <= Ni; j++) { 
     grad[i] -= 1.0 / j;
    }
    grad[imin] -= grad[i];
  }

}
*/

double FuncionEntropiaNatural::f(double pars[]) {
  int i;
  double S = 0, N = 0;

  for (i = 0; i < n_pars; i++) {
    N += pars[i];
    S -= lgamma(pars[i] + 1);
  }
  S += lgamma (N + 1) - N * log ((double) n_pars);

  if (isinf(S)) {
    S = 0;
  }
  
  return S;
}


void FuncionEntropiaNatural::df(double pars[], double grad[]){
  int i, j;
  double dS_cte = 0, N = 0;
  
  for (i = 0; i < n_pars; i++) {
    N += pars[i];
  }

  dS_cte = -log ((double) n_pars);
  for (i = 1; i <= N; i++) {
    dS_cte += 1.0 / i;
  }
  
  for (i = 0; i < n_pars; i++) {
    grad[i] = dS_cte;
    for (j = 1; j <= pars[i]; j++) {
     grad[i] -= 1.0 / j;
    }
  }
}


