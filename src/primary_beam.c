#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mathconst.h"
#include "primary_beam.h"


/*--------------------------------------------------------------------
 * Initialize the parameters of the primary beam.
 *--------------------------------------------------------------------*/

void init_beam(struct pbeam *beam)
{
  beam->type = GAUSS;
  beam->fwhm = 44.0*RPARCM;   /* radians */
  beam->freq = 30.0;          /* GHz */
  beam->cutoff = 90.0*RPARCM; /* radians */
}

/*--------------------------------------------------------------------
 * Return the response of the primary beam at angle arc radian from 
 * pointing direction at frequency freq GHz.
 *--------------------------------------------------------------------*/

double primary_beam(double arc, double freq, struct pbeam *beam)
{
  const double c = 4.0*LN2; /* 4 ln(2) */
  double a, r;

  if (arc > beam->cutoff) {
    return 0.0;
  }
  switch (beam->type) {
  case GAUSS:
    /* Gaussian approximation */
    a = (beam->fwhm)*((beam->freq)/freq); /* FWHM in radians */
    r = arc/a;
    return exp(-c*r*r);
    break;
  case CBI:
    return cbi_beam(arc, freq);
    break;
  case DESIGN:
    return cbi_beam_design(arc, freq);
    break;
  case CBI2:
    return cbi2_beam(arc, freq);
    break;
  default:
    fprintf(stderr, "Internal error: bad beam type\n");
    exit(EXIT_FAILURE);
    break;
  }
}
