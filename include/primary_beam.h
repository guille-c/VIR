#include "cbi_beam.h"

#ifdef __cplusplus
extern "C" {
#endif  

/* Structure to hold primary beam parameters */

enum BEAMTYPE {GAUSS = 0, CBI = 1, DESIGN = 2, CBI2 = 3};
struct pbeam {
  enum BEAMTYPE type;    /* Functional form */
  double fwhm;           /* Full-width at half maximum (radian) (GAUSS) */
  double freq;           /* Frequency at which this FWHM applies (GHz) */
  double cutoff;         /* Maximum radius in radian */
};

void init_beam(struct pbeam *beam);
double primary_beam(double arc, double freq, struct pbeam *beam);

#ifdef __cplusplus
}
#endif  

