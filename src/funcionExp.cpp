#include "funcionExp.h"

double FuncionExp::f(double pars[]) {
  double ret = 0;
  for (int i = 0; i < n_pars; i++) {
    ret += 0.4 * pow (pars[i] - 0.4, 2) - 0.08 * pow (pars[i], 4);
  }
  //cout << "ret = " << ret << ", f = " << exp(ret) << "\n";
  return exp (ret / n_pars);
}

void FuncionExp::df(double pars[], double grad[]) {
  for (int i = 0; i < n_pars; i++) {
    grad[i] = f(pars) * (0.8 * (pars[i] - 0.4) - 0.32 * pow (pars[i], 3)) / n_pars;
    //cout << "p[" << i << "] = " << pars[i] << ", grad[" << i << "] = " << grad[i] 
    //<< ", f = " << f(pars) << "\n";
  }
}
