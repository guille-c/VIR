#ifndef DISTRIBUCIONCONDICIONALCBI
#define DISTRIBUCIONCONDICIONALCBI

#include <math.h>
#include <string>
#include <fstream.h>
#include "distribucionCondicional.h"
#include "mockcbiRoutines.h"

#define SQR(x)  ((x)*(x))   // Funcion para calcular cuadrados.

class DistribucionCondicionalCBI: public DistribucionCondicional {
 protected:
  char *nombreImagen, **nombresVis;
  int nVis;
  double ruido, ref_freq, fg_scale;
  struct image *fg_image, *cmb_image;
  struct uvf_header **header_obs;  // Arreglo con los encabezados.
  struct uvf_sample **samples_obs; // Arreglo con las visibilidades observadas.
  struct uvf_sample **samples_mod; // Arreglo con las visibilidades modeladas.
  struct pbeam beam;               // Beam del CBI
  double ***atten;

 public:
  DistribucionCondicionalCBI (char *nombreImagen, char **nombresVis, int nVis);
  double getQuanta() {return ruido;}

  virtual void guardarInfo (string s);

  virtual ~DistribucionCondicionalCBI();
};

#endif
