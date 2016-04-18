#ifndef DISTRIBUCIONPBB
#define DISTRIBUCIONPBB

#include <string>
using namespace std;

class DistribucionPbb {
 protected:
  int n_pars;
 public:
  DistribucionPbb(int n_pars) {this->n_pars = n_pars;}
  virtual double f(double *pars) = 0;
  virtual void generar(double *pars) = 0;
  virtual void guardarInfo(string s) {}
  int getNPars() {return n_pars;}
};

#endif
