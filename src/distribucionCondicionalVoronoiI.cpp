#include "distribucionCondicionalVoronoiI.h"

DistribucionCondicionalVoronoiI::DistribucionCondicionalVoronoiI(DistribucionCondicional *dc, Inicializador *ini) {
  if (dc->getNPars() % 3 != 0) {
    cerr << "ERROR al intentar inicializar DistribucionCondicionalVoronoiI\n";
    cerr << "      dc no tiene 3*n parametros. \n";
    exit (1);
  }
  this->dc = dc;
  n_pars = dc->getNPars() / 3;
  pos = new double [n_pars * 2];
  
  double *pos_aux = new double [3 * this->n_pars];
  ini->inicializar (pos_aux, this->n_pars * 3);
  for (int i = 0; i < this->n_pars; i++) {
    pos [2 * i] = pos_aux [3 * i];
    pos [2 * i + 1] = pos_aux [3 * i + 1];
  }
  delete [] pos_aux;
}

double DistribucionCondicionalVoronoiI::generar (double *p, int j) {
  
  cout << "generando 1 " << j << "\n";
  double *p_aux = new double [3 * this->n_pars];
  for (int i = 0; i < this->n_pars; i++) {
    p_aux [3 * i] = pos [2*i];
    p_aux [3 * i + 1] = pos [2*i + 1];
    p_aux [3 * i + 2] = p [i];
  }
  cout << "generando 2 " << j << "\n";
  double ret = dc->generar (p_aux, 3*j + 2);
  cout << "generado " << j << "\n";
  delete [] p_aux;
  return ret;
}

double DistribucionCondicionalVoronoiI::evaluar (double *p) {
  double *p_aux = new double [3 * this->n_pars];
  for (int i = 0; i < this->n_pars; i++) {
    p_aux [3 * i] = pos [2*i];
    p_aux [3 * i + 1] = pos [2*i + 1];
    p_aux [3 * i + 2] = p [i];
  }

  double ret = dc->evaluar (p_aux);
  delete [] p_aux;
  return ret;
}

double DistribucionCondicionalVoronoiI::evaluar (double *p, double x, int j) {
  cout << "evaluando 1 " << j << "\n";
  double *p_aux = new double [3 * this->n_pars];
  for (int i = 0; i < this->n_pars; i++) {
    p_aux [3 * i] = pos [2*i];
    p_aux [3 * i + 1] = pos [2*i + 1];
    p_aux [3 * i + 2] = p [i];
  }

  cout << "evaluando 2 " << j << "\n";
  double ret = dc->evaluar (p_aux, x, 3*j + 2);
  cout << "evaluado " << j << "\n";
  delete [] p_aux;
  return ret;
}
