#include "funcionReconstructor.h"

FuncionReconstructor::FuncionReconstructor (Reconstructor *rec) {
  this->n_pars = rec->getNPars();
  this->rec = rec; 
}

double FuncionReconstructor::f(double pars[]) {
  double ret;
  for (int i = 0; i < n_pars; i++) {
    //cout << i << ") Cambiando " << rec->getPar(i) << " por " << pars[i] << "\n";
    rec->setPar (i, pars[i]);
  }
  ret = rec->run();

  for (int i = 0; i < n_pars; i++) {
    //cout << i << ") Cambiando " << pars[i] << " por " << rec->getPar (i) << "\n";
    pars[i] = rec->getPar (i);
  }
  return ret;
}

void FuncionReconstructor::df(double pars[], double grad[]){
  // No Se :P
}

void FuncionReconstructor::guardarInfo (string info) {
  string s = "funcRec_";
  rec->guardarInfo (s + info);
}

void FuncionReconstructor::setNPars (int n_pars) {
  double *p_aux = new double [n_pars];
  int n = n_pars;
  if (n > this->rec->getNPars()) {
    n = this->rec->getNPars();
  }

  this->n_pars = n_pars;
  for (int i = 0; i < n; i++) {
    p_aux[i] = this->rec->getPar (i);
  }
  this->rec->setPars(p_aux);
  this->rec->setNPars (n_pars);
}

FuncionReconstructor::~FuncionReconstructor() {
  if (rec != NULL) {
    delete rec;
  }
  rec = NULL;
}
