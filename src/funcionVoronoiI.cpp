#include "funcionVoronoiI.h"

FuncionVoronoiI::FuncionVoronoiI (FuncionVoronoi *fV, Inicializador *ini) {
  this->fV = fV;
  this->n_pars = fV->getNPars() / 3;

  pos = new double [2 * fV->getNPars()];

  double *pos_aux = new double [3 * this->n_pars];
  ini->inicializar (pos_aux, this->n_pars * 3, 0.0);
  for (int i = 0; i < this->n_pars; i++) {
    pos [2 * i] = pos_aux [3 * i];
    pos [2 * i + 1] = pos_aux [3 * i + 1];
  }
  delete [] pos_aux;
}

FuncionVoronoiI::FuncionVoronoiI (FuncionVoronoi *fV, double *pos) {
  this->fV = fV;
  this->n_pars = fV->getNPars() / 3;

  this->pos = new double [2 * fV->getNPars()];

  for (int i = 0; i < this->n_pars * 2; i++) {
    this->pos [i] = pos [i];
  }
}

double FuncionVoronoiI::f (double *pars) {
  double ret, *parsF;
  parsF = new double [this->n_pars * 3];
  for (int i = 0; i < this->n_pars; i++) {
    parsF [3 * i] = pos[2 * i];
    parsF [3 * i + 1] = pos[2 * i + 1];
    parsF [3 * i + 2] = pars[i];
  }
  ret = fV->f(parsF);
  delete [] parsF;
  return ret;
}

void FuncionVoronoiI::df (double *pars, double * grad) {
  double ret, *parsF, *gradF;

  parsF = new double [this->n_pars * 3];
  gradF = new double [this->n_pars * 3];
  for (int i = 0; i < this->n_pars; i++) {
    parsF [3 * i] = pos[2 * i];
    parsF [3 * i + 1] = pos[2 * i + 1];
    parsF [3 * i + 2] = pars[i];
  }
  fV->df (parsF, gradF);

  for (int i = 0; i < this->n_pars; i++) {
   grad[i] = gradF [3 * i + 2];
  }

  delete [] parsF;
  delete [] gradF;
}

void FuncionVoronoiI::guardarInfo (string info) {
  fV->guardarInfo (info);
}

void FuncionVoronoiI::setNPars (int n_pars) {
  this->n_pars = n_pars;
  fV->setNPars (n_pars * 3);
  if (pos != NULL) {
    delete [] pos;
    pos = new double [n_pars * 2];
  }
}

FuncionVoronoiI::~FuncionVoronoiI() {
  if (pos != NULL) {
    delete [] pos;
    pos = NULL;
  }
}
