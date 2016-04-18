#ifndef VIRPP_FUNCIONCBI
#define VIRPP_FUNCIONCBI

#include "funcion.h"
#include "uvsubs.h"
//#include "image_routines.h"
//#include "primary_beam.h"
#include "mockcbiRoutines.h"

#define SQR(x)  ((x)*(x))   // Funcion para calcular cuadrados.

class FuncionCBI: public virtual Funcion {
 protected:
  char *nombreImagen, **nombresVis;
  int nVis;
  double ruido, ref_freq, fg_scale;
  struct image *fg_image, *cmb_image;
  struct uvf_header **header_obs;  // Arreglo con los encabezados.
  struct uvf_sample **samples_obs; // Arreglo con las visibilidades observadas.
  struct uvf_sample **samples_mod; // Arreglo con las visibilidades modeladas.
  struct pbeam beam;               // Beam del CBI
  struct image *noise_image;

  double ***atten;
 public:
  FuncionCBI (): Funcion (){}
  FuncionCBI (char * nombreImagen, char **nombresVis, 
	      int nVis, int n_pars = 0, double beamSize = -1);
  virtual double getRuido () {return ruido;}
  virtual void setRuido (double r) {ruido = r;}
  virtual int getNDat () {return 2*nVis * header_obs[0]->nif * header_obs[0]->nsamp;}
  virtual void resizeImage (int nx, int ny, int cx, int cy);
  double cross_correlate_factor ();
  void guardarResiduos (string info);
  void calcularNoiseImage();

  virtual ~FuncionCBI ();

  FuncionCBI& operator= (const FuncionCBI&);
};

#endif
