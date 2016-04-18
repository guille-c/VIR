#ifndef DISTRIBUCIONCONDICIONAL
#define DISTRIBUCIONCONDICIONAL

#include <string>
#include <math.h>
#include <iostream.h>
using namespace std;

#define DOUBLE_MAX 1.79769E+308
#define DOUBLE_MIN -1.79769E+308

class DistribucionCondicional {
 protected:
  int n_pars;
  
 public:
  virtual double generar (double *, int i) = 0;
  virtual double evaluar (double *) = 0;
  virtual double evaluar (double *, double x, int i) = 0;
  virtual int getNPars() {return n_pars;}
  
  virtual void busquedaCotas (double *p, int i, double *a, double *b);
  virtual void guardarInfo (string s) {}
};

#endif
