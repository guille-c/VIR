#include "reconstructorIndividual.h"

ReconstructorIndividual::ReconstructorIndividual (Funcion *f, double min, double max, double N, Inicializador *ini)
  : Reconstructor (f->getNPars(), ini) {
  this->f = f;
  this->N = new double[n_pars];
  this->max = new double[n_pars];
  this->min = new double[n_pars];
  for (int i = 0; i < n_pars; i++) {
    this->max [i] = max;
    this->min [i] = min;
    this->N[i] = N;
  }
}

ReconstructorIndividual::ReconstructorIndividual (Funcion *f, double *min, double *max, double *N, Inicializador *ini)
  : Reconstructor (f->getNPars(), ini) {
  this->f = f;
  this->N = new double[n_pars];
  this->max = new double[n_pars];
  this->min = new double[n_pars];
  for (int i = 0; i < n_pars; i++) {
    this->max [i] = max[i];
    this->min [i] = min[i];
    this->N[i] = N[i];
  }
}

double ReconstructorIndividual::run () {
  double dx, fmax, f, pmax;
  
  fmax = this->f->f(p);
  for (int i = 0; i < n_pars; i++) {
    cout << "Buscando " << i << " \n";
    dx = (max[i] - min[i]) / N[i];
    //fmax = -1e300;
    for (double x = min[i]; x <= max[i]; x += dx) {
      //cout << "x = " << x << "\n";
      p[i] = x;
      f = this->f->f(p);
      if (f > fmax) {
	fmax = f;
	pmax = x;
      }
    }
    p[i] = pmax;
    cout << "p[" << i << "] = " << p[i] << "\n";
    this->f->f(p);
    this->f->guardarInfo(i2s(i, 3));
  }
  return fmax;
}

ReconstructorIndividual::~ReconstructorIndividual() {
  delete [] max;
  delete [] min;
}
