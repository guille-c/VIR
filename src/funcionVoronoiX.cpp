#include "funcionVoronoiX.h"

FuncionVoronoiX::FuncionVoronoiX (FuncionVoronoi *fV, Inicializador *ini) {
  this->fV = fV;
  this->n_pars = fV->getNPars() / 3 * 2;

  I = new double [fV->getNPars() / 3];

  double *pos_aux = new double [fV->getNPars()];
  ini->inicializar (pos_aux, this->n_pars * 3 / 2, 0.0);
  for (int i = 0; i < fV->getNPars() / 3; i++) {
    I [i] = pos_aux [3 * i + 2];
  }
  delete [] pos_aux;
}

FuncionVoronoiX::FuncionVoronoiX (FuncionVoronoi *fV, double *I) {
  this->fV = fV;
  this->n_pars = fV->getNPars() / 3 * 2;

  this->I = new double [fV->getNPars() / 3];

  for (int i = 0; i < fV->getNPars() / 3; i++) {
    this->I [i] = I[i];
  }
}

double FuncionVoronoiX::f (double *pars) {
  double ret, *parsF;
  parsF = new double [this->n_pars * 3 / 2];
  for (int i = 0; i < this->n_pars / 2; i++) {
    parsF [3 * i] = pars[2 * i];
    parsF [3 * i + 1] = pars[2 * i + 1];
    parsF [3 * i + 2] = I[i];
  }
  ret = fV->f(parsF);
  delete [] parsF;
  return ret;
}

void FuncionVoronoiX::df (double *pars, double *grad) {
  double ret, *parsF, *gradF;

  parsF = new double [this->n_pars * 3 / 2];
  gradF = new double [this->n_pars * 3 / 2];
  for (int i = 0; i < this->n_pars / 2; i++) {
    parsF [3 * i] = pars[2 * i];
    parsF [3 * i + 1] = pars[2 * i + 1];
    parsF [3 * i + 2] = I[i];
  }
  fV->df (parsF, gradF);

  for (int i = 0; i < this->n_pars / 2; i++) {
   grad[2 * i] = gradF [3 * i];
   grad[2 * i + 1] = gradF [3 * i + 1];
  }

  delete [] parsF;
  delete [] gradF;
}

void FuncionVoronoiX::guardarInfo (string info) {
  fV->guardarInfo (info);
}

void FuncionVoronoiX::setNPars (int n_pars) {
  this->n_pars = n_pars;
  fV->setNPars (n_pars * 3 / 2);
  if (I != NULL) {
    delete [] I;
    I = new double [n_pars / 2];
  }
}

FuncionVoronoiX::~FuncionVoronoiX() {
  if (I != NULL) {
    delete [] I;
    I = NULL;
  }
}
