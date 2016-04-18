#include "distribucionNormal.h"

DistribucionNormal::DistribucionNormal (int n_pars, double sigma, int seed, int n)
  : DistribucionPbb (n_pars){
  this->n = n;
  this->media = new double [n_pars];
  for (int i = 0; i < n_pars; i++) {
    this->media[i] = 0;
  }  
  this->sigma = sigma;
  this->add = sqrt(3 * n);
  this->fac = 2 * add / ( (double ) n * RAND_MAX);
  srand (seed);
}

DistribucionNormal::DistribucionNormal (int n_pars, double *media, double sigma, int seed, int n)
  : DistribucionPbb (n_pars){
  this->n = n;
  this->media = new double [n_pars];
  for (int i = 0; i < n_pars; i++) {
    this->media[i] = media[i];
  }
  this->sigma = sigma;
  this->add = sqrt(3 * n);
  this->fac = 2 * add / ( (double ) n * RAND_MAX);
  srand (seed);
}

double DistribucionNormal::f(double pars[]) {
  double suma = 0;

  for (int i = 0; i < n; i++) {
    suma += -(pars[i] - media[i]) * (pars[i] - media[i]) / (2 * sigma * sigma);
  }
  return exp (suma) / (sigma * sqrt (2 * PI));
}

void DistribucionNormal::generar(double pars[]) {
  double sum;
  int i, j;
  
  for (j = 0; j < n_pars; j++) {
    sum = 0;
    for (i = 0; i < n; i++)
      sum += rand();
    pars[j] = (fac * sum - add) * sigma + media[j];
  }
}

void DistribucionNormal::setMedia (double *media) {
  for (int i = 0; i < n_pars; i++) {
    this->media[i] = media[i];
  }
}

void DistribucionNormal::setSigma (double sigma) {
  this->sigma = sigma;
}

DistribucionNormal::~DistribucionNormal () {
  delete [] media;
}
