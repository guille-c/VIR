#include "distribucionCondicionalBVCBI_I.h"

DistribucionCondicionalVoronoiI::DistribucionCondicionalVoronoiI(DistribucionCondicional *dc, Inicializador *ini) {
  if (dc->getNPars() % 3 != 0) {
    cerr << "ERROR al intentar inicializar DistribucionCondicionalVoronoiI\n";
    cerr << "      dc no tiene 3*n parametros. \n";
    exit (1);
  }
  n_pars = f->getNPars();
  this->N = N;
}
