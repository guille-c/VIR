#ifndef DISTRIBUCIONCONDGAUSSIANA
#define DISTRIBUCIONCONDGAUSSIANA

#include <math.h>
#include "distribucionCondicional.h"
#include "distribucionNormal.h"

class DistribucionCondGaussiana: public DistribucionCondicional {
 protected:
  double *medias;
  double *sigmas;
  DistribucionNormal *distN;
 public:
  DistribucionCondGaussiana (int n_pars, double *medias, double *sigmas);
  virtual double generar (double *, int i);
  virtual double evaluar (double *);
  virtual double evaluar (double *, double x, int i);

  virtual ~DistribucionCondGaussiana();
};

#endif
