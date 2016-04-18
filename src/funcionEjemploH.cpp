#include "funcionEjemploH.h"

double FuncionEjemploH::f(double pars[]) {
  double ret = 0;
  for (int i = 0; i < n_pars; i++) {
    ret += 0.4 * pow (pars[i] - 0.4, 2) - 0.08 * pow (pars[i], 4);
  }
  //return exp (ret);
  return -ret;
}

void FuncionEjemploH::df(double pars[], double grad[]) {
  //double func = f(pars);
  for (int i = 0; i < n_pars; i++) {
    //grad[i] = func * (0.8 * (pars[i] - 0.4) - 0.32 * pow (pars[i], 3));
    grad[i] = -(0.8 * (pars[i] - 0.4) - 0.32 * pow (pars[i], 3));
  }
}
