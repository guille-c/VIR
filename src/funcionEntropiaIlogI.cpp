#include "funcionEntropiaIlogI.h"

double FuncionEntropiaIlogI::f(double pars[]) {
  int i;
  double S = 0;

  if (lambda == 0) return 0;

  for (i = 0; i < n_pars; i++) {
    S -= pars[i] * log (pars[i] / M);
  }

  return lambda * S;
}


void FuncionEntropiaIlogI::df(double pars[], double grad[]){
  int i;
  
  if (lambda == 0) {
    for (i = 0; i < n_pars; i++) {
      grad[i] = 0;
    }
    return;
  }

  for (i = 0; i < n_pars; i++) {
    grad[i] = -lambda * (1 + log(pars[i] / M));
  }

}


