#include "funcionDividida.h"

FuncionDividida::FuncionDividida (Funcion *f) {
  funcion = f; 
  this->n_pars = f->getNPars();
}

double FuncionDividida::f(double pars[]) {
  return 1.0 / funcion->f(pars);
}

void FuncionDividida::df(double pars[], double grad[]) {
  funcion->df(pars, grad);
  for (int i = 0; i < n_pars; i++) {
    grad[i] *= -1 / pow (f(pars), 2);
  }
}

void FuncionDividida::setNPars (int n_pars) {
  this->n_pars = n_pars;
  funcion->setNPars (n_pars);
}

void FuncionDividida::guardarInfo (string info) {
  funcion->guardarInfo(info);
}

FuncionDividida::~FuncionDividida() {
  if (funcion != NULL) {
    delete funcion;
    funcion = NULL;
  }
}
