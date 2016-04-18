#ifndef VIRPP_FUNCION
#define VIRPP_FUNCION

#include <iostream.h>
#include <fstream.h>
#include <string>
using namespace std;

class Funcion {
 protected:
  int n_pars;
 public:
  Funcion () {n_pars = 0;}
  Funcion (int n) {this->n_pars = n;}
  virtual double f(double pars[]) = 0;
  virtual void df(double pars[], double grad[]) = 0;
  int getNPars() {return n_pars;}
  void setNPars (int n_pars) {this->n_pars = n_pars;}
  virtual void guardarInfo (string info){}

  virtual ~Funcion() {}
};

#endif
