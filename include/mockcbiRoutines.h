#ifndef MOCKCBI_ROUTINES
#define MOCKCBI_ROUTINES

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "complex.h"        /* Complex arithmetic */
#include "fitsio.h"         /* CFITSIO library */
#include "fourn.h"          /* FFT routine fourn() */
#include "image_routines.h" /* For reading FITS images */
#include "mathconst.h"      /* Constants */
#include "parse_sexagesimal.h" /* defines parse_sexagesimal() */
#include "primary_beam.h"   /* define primary_beam(), init_beam() */
#include "ra_dec_string.h"
#include "random.h"         /* defines gasdev() */
#include "slalib.h"         /* Patrick Wallace's SLA library */
#include "uvsubs.h"         /* UVFITS routines */

#ifdef __cplusplus
extern "C" {
#endif  

struct leakage {
  double real;
  double imag;
  int defined;
};
struct leakage_table {
  int do_leakage;
  int n_rx;
  int n_ch;
  struct leakage **leakage;
};

/* Status returns */

typedef enum {SUCCESS, FAILURE} Status;

double ***atenuacion(struct image *fg_image, struct uvf_header **header_obs,
		     int n_archivos, struct pbeam beam);
void mockcbi_sub(struct image *cmb_image, 
		 struct image *fg_image,
		 struct uvf_header *header, 
		 struct uvf_sample *samples,
		 struct uvf_sample *samples_salida,
		 struct pbeam beam);
void copiar_uvf_samples(struct uvf_sample *sin, struct uvf_sample *sout, int n);
void multiplicar_vis(struct uvf_header *header,
		     struct uvf_sample *samples, 
		     double factor);
Status do_clear(char *command, struct uvf_header *header,
		struct uvf_sample *vis);
Status do_read_uvdata(char *file, struct uvf_header **header,
		      struct uvf_sample **samples);
Status do_write(char *command, char *file, struct uvf_header *header,
		struct uvf_sample *samples, char *infile);
Status do_add_cmb(char *command, double scalefactor,
		  double ra_offset, struct uvf_header *header,
		  struct uvf_sample *samples,
		  struct image *cmb_image,
		  struct image *q_image,
		  struct image *u_image,
		  struct image *fg_image, double fg_spind,
		  struct pbeam *beam,
		  struct leakage_table *q);
void apply_beam(float *image, struct image *sky,
		struct image *fg, double fg_spind,
		double x0, double y0, double freq,
		struct pbeam *beam);
Status getvis(double u, double v, float *image, long nx, long ny,
	      float *visr, float *visi);
void phase_rotate(float *image, long nx, long ny,
		  double x, double y);
void direccos(double ra, double dec, double ra0, double dec0,
	      double *l, double *m);
int check_fg(struct image *img1, struct image *img2);

void scale (struct uvf_header *header,             /* scale visibilities  */
	    struct uvf_sample *samples, float ascale) ;
void rotar_vis(struct uvf_header *header,
	       struct uvf_sample *samples);
double Irms (struct uvf_header **header,
	     struct uvf_sample **samples, int n_archivos);
void residual_vis(struct uvf_header *header,
			 struct uvf_sample *samples,
			 struct uvf_sample *samples_mod,  
			 struct uvf_sample *samples_res);
void dirtyMap (struct image *im, struct uvf_header **header, struct uvf_sample **samples, 
	       int nVis, double ***atten);

#ifdef __cplusplus
}
#endif  

#endif
