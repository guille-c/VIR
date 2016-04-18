#include "distribucionCondGaussiana.h"

DistribucionCondGaussiana::DistribucionCondGaussiana (int n_pars, double *medias, double *sigmas){
  this->n_pars = n_pars;
  this->medias = new double [n_pars];
  this->sigmas = new double [n_pars];
  for (int i = 0; i < n_pars; i++) {
    this->medias[i] = medias[i];
    this->sigmas[i] = sigmas[i];
    cout << i <<") " << this->medias[i] << ", " << this->sigmas[i] << "\n";
  }
  double *media = new double [1];
  media[0] = medias[0];
  distN = new DistribucionNormal (1, media, sigmas[0]);
  delete [] media;
}

double DistribucionCondGaussiana::generar (double *p, int i) {
  double *media = new double[1];
  double *pars = new double[1];
  double ret;
  media[0] = medias[i];
  distN->setMedia (media);
  distN->setSigma(sigmas[i]);
  distN->generar(pars);
  ret = pars[0];
  //cout << ret << ", " << media[0] << ", " << sigmas[i] << "\n";
  delete [] media;
  delete [] pars;

  return ret;
}
 
double DistribucionCondGaussiana::evaluar (double *p) {
  double ret = 1;
  
  for (int i = 0; i < n_pars; i++) {
    evaluar (p, p[i], i);
  }
  return ret;
}

double DistribucionCondGaussiana::evaluar (double *p, double x, int i) {
  double *media = new double[1];
  double *pars = new double[1];
  double ret = 1;
  double x_old = p[i];
  p[i] = x;

  for (int j = 0; j < n_pars; j++) {
    media[0] = medias[j];
    distN->setMedia (media);
    distN->setSigma (sigmas[j]);
    pars[0] = p[j];
    ret *= distN->f(pars);
  }
  p[i] = x_old;
  delete [] media;
  delete [] pars;

  return ret;
}

DistribucionCondGaussiana::~DistribucionCondGaussiana() {
  delete [] medias;
  delete [] sigmas;
  delete distN;
}
