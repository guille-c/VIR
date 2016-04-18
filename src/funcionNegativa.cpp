#include "funcionNegativa.h"

FuncionNegativa::FuncionNegativa (Funcion *f) {
  funcion = f; 
  this->n_pars = f->getNPars();
}

double FuncionNegativa::f(double pars[]) {
  return -funcion->f(pars);
}

void FuncionNegativa::df(double pars[], double grad[]) {
  funcion->df(pars, grad);
  for (int i = 0; i < n_pars; i++) {
    grad[i] *= -1;
  }
}

void FuncionNegativa::setNPars (int n_pars) {
  this->n_pars = n_pars;
  funcion->setNPars (n_pars);
}

void FuncionNegativa::guardarInfo (string info) {
  funcion->guardarInfo(info);
}

FuncionNegativa::~FuncionNegativa() {
  if (funcion != NULL) {
    delete funcion;
    funcion = NULL;
  }
}
