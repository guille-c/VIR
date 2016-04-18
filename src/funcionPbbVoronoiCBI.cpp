#include "funcionPbbVoronoiCBI.h"

FuncionPbbVoronoiCBI::FuncionPbbVoronoiCBI (char * nombreImagen, char **nombresVis, 
					    int nVis, int n,
					    int entropia, double cuantaSize,
					    double expanded_lx, double expanded_ly){

  fBV = new FuncionBayesVoronoiCBI (nombreImagen, nombresVis, nVis,  n,
				    entropia, cuantaSize, expanded_lx, expanded_ly);
  this->n_pars = fBV->getNPars();
  this->malla = fBV->getMalla();
  double *pars = new double[n_pars];

  for (int i = 0; i < n_pars / 3; i++) {
    pars[3 * i]  = rand()/(RAND_MAX + 1.0);
    pars[3 * i + 1]  = rand()/(RAND_MAX + 1.0);
    pars[3 * i + 2]  = 0;
  }
  this->c = fabs (fBV->f(pars));
  delete [] pars;
}

FuncionPbbVoronoiCBI::FuncionPbbVoronoiCBI (char *nombre_archivo) {

  fBV = new FuncionBayesVoronoiCBI (nombre_archivo);
  this->n_pars = fBV->getNPars();
  this->malla = fBV->getMalla();
  double *pars = new double[n_pars];

  for (int i = 0; i < n_pars / 3; i++) {
    pars[3 * i]  = rand()/(RAND_MAX + 1.0);
    pars[3 * i + 1]  = rand()/(RAND_MAX + 1.0);
    pars[3 * i + 2]  = 0;
  }
  this->c = fabs (fBV->f(pars));
  delete [] pars;
}

double FuncionPbbVoronoiCBI::f (double *pars) {

  for (int i = 0; i < n_pars / 3; i++) {
    if (pars[3*i] < 0 || pars[3*i] >= 1 ||
	pars[3*i + 1] < 0 || pars[3*i + 1] >= 1) {
      return 1e-300;
    }
  }
  double ret = exp(-(fBV->f(pars)/c));

  if (ret > 1) {
    cerr << "ERROR: en f, valor " << ret << " > 1\n";
    exit (1);
  }

  return ret;
}

void FuncionPbbVoronoiCBI::df (double *pars, double *grad){
  double func = f(pars);

  fBV->df(pars, grad);
  for (int i = 0; i < n_pars; i++) {
    grad[i] *= -func/c;
  }
}

void FuncionPbbVoronoiCBI::guardarInfo (char *info){
  // Para este caso usaremos info como el sufijo de los archivos a
  // guargar.
  
  cout << "Guardando info " << info << "\n";
  fBV->guardarInfo(info);
}

void FuncionPbbVoronoiCBI::guardarInfo (string info){
  // Para este caso usaremos info como el sufijo de los archivos a
  // guargar.
  
  cout << "Guardando info " << info << "\n";
  fBV->guardarInfo(info);
}


FuncionPbbVoronoiCBI::~FuncionPbbVoronoiCBI() {
  delete fBV;
}
