#ifndef DISTRIBUCIONNORMAL
#define DISTRIBUCIONNORMAL

#include "distribucionPbb.h"
#include <math.h>
#include <cstdlib>

#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751
#endif

class DistribucionNormal: public DistribucionPbb {
 protected:
  int n;
  double *media, sigma;
  double add, fac;  // variables para calcular la distribucion.
 public:
  DistribucionNormal (int n_pars, double sigma, int seed = 0, int n = 10);
  DistribucionNormal (int n_pars, double *media, double sigma, int seed = 0, int n = 10);
  double f (double *pars);
  void generar (double *pars);
  void setMedia (double *media);
  void setSigma (double sigma);

  ~DistribucionNormal ();
};

#endif
