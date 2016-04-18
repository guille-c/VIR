#ifndef DISTRIBUCIONCONDICIONALBAYESVORONOICBI
#define DISTRIBUCIONCONDICIONALBAYESVORONOICBI

#include <math.h>
#include <string>
#include "distribucionCondicional.h"
#include "funcionBayesVoronoiCBI.h"

class DistribucionCondicionalBVCBI: public DistribucionCondicional {
 protected:
  FuncionBayesVoronoiCBI *f;
  int N;
 public:
  DistribucionCondicionalBVCBI (char *nombreArchivoInput, int N);
  DistribucionCondicionalBVCBI(char * nombreImagen, char **nombresVis, 
			       int nVis, int n,
			       int entropia, double cuantaSize,
			       double expanded_lx, double expanded_ly, 
			       int N);
  virtual double generar (double *, int i);
  double generarIntensidad (double *, int i);
  double generarPosicion (double *, int i);
  virtual double evaluar (double *);
  virtual double evaluar (double *, double x, int i);
  FuncionBayesVoronoiCBI* getFunc() {return f;}
  double getQuanta() {return f->getRuido();}

  virtual void guardarInfo (string s) {f->guardarInfo(s);}
  void busquedaCotas (double *p, int i, double *a, double *b);

  virtual ~DistribucionCondicionalBVCBI();
};

#endif
