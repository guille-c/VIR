/* uvsubs.h */

#ifndef UVSUBS
#define UVSUBS

#include "fitsio.h"

#ifdef __cplusplus
extern "C" {
#endif  
#define MAX_DIM 100    /* Maximum number of axes in a UV FITS file */
#define MAX_PAR 20     /* Maximum number of "random parameters" 
			  in a UV FITS file */
#define MAX_ANTENNAS 100 /* Maximum number of antennas in an antenna table */
#define MAX_DAT 120
#define MAX_STOKES 4   /* Maximum number of pixels on STOKES axis */

/* Structure describing an antenna */

struct ant {
  char name[32];                    /* Antenna name */
  double x;                         /* x, y, z coordinates */
  double y;
  double z;
};

/* Structure for holding a file header; this contains all
 * the information about the dataset except for the number of
 * samples and the samples themselves. */

struct uvf_header {
  double equinox;                   /* Equinox of coordinates, usually 2000 */
  double obsra;                     /* Antenna pointing RA [deg] */
  double obsdec;                    /* Antenna pointing Dec [deg] */
  double ra;                        /* Phase center RA [deg] */
  double dec;                       /* Phase center Dec [deg] */
  double freq;                      /* Reference frequency [Hz] */
  double bw;                        /* Channel bandwidth [Hz] */
  double start_jd;                  /* Start time */
  double end_jd;                    /* End time */
  double pscal[MAX_PAR];            /* Scale factors for parameters */
  double pzero[MAX_PAR];            /* Offsets for parameters */
  double iffreq[MAX_DAT];           /* Frequency of each IF */
  int nstokes;                      /* Number of polarizations */
  int stokes[MAX_STOKES];           /* Stokes parameters */
  int delta_stokes;                 /* Axis increment in Stokes */
  int nif;                          /* Number of IFs (frequency channels) */
  int nchan;                        /* Number of complex visibilities in each
				       group */
  int nsamp;                        /* Number of samples */
  int pcount;                       /* Number of parameters */
  int index_u;                      /* Parameter number for u */
  int index_v;                      /* Parameter number for v */
  int index_w;                      /* Parameter number for w */
  int index_baseline;               /* Parameter number for baseline */
  int index_date1;                  /* Parameter number for date (1) */
  int index_date2;                  /* Parameter number for date (2) */
  int index_inttim;                 /* Parameter number for integration time */
  char object[FLEN_VALUE];          /* Object (source) name */
  char telescop[FLEN_VALUE];        /* Telescope name */
  char instrume[FLEN_VALUE];        /* Instrument name */
  char bunit[FLEN_VALUE];           /* Units of data */
  char radecsys[FLEN_VALUE];        /* Coordinate system, usually "FK5" */
  char observer[FLEN_VALUE];        /* Observer name or code */
  
  int n_antennas;                   /* Number of antennas */
  struct ant antenna[MAX_ANTENNAS]; /* Info on each antenna */
};

/* Structure for holding a single visibility sample.  One sample
 * is characterized by baseline number, u,v,w, date and nif
 * complex number with weights (3*nif numbers in all). */

struct uvf_sample {
  double date;                      /* Modified Julian Date */
  double u;                         /* Baseline u component [sec] */
  double v;                         /* Baseline v component [sec] */
  double w;                         /* Baseline w component [sec] */
  double inttim;                    /* Integration time */
  float rdata[MAX_DAT];             /* Complex visibilities and weights */
  int baseline;                     /* Baseline number */
};

/* Routines defined in this file */

int get_header(fitsfile *fptr, struct uvf_header *header);
int get_sample(fitsfile *fptr, struct uvf_header *header,
	       int number, struct uvf_sample *sample);
int get_antable(fitsfile *fptr, int hdu, struct uvf_header *header);
int get_fqtable(fitsfile *fptr, int hdu, struct uvf_header *header);
int put_header(fitsfile *fptr, struct uvf_header *header, int samples,
	       double refdat, char *origin);
int put_sample(fitsfile *fptr, struct uvf_header *header,
	       int group, double refdat, struct uvf_sample *sample);
int comp_headers(struct uvf_header *h1, struct uvf_header *h2,
		 int f1, int f2);
int copy_history(fitsfile *infptr, fitsfile *outfptr, char *infile);
char *stokes_label(int code);

#ifdef __cplusplus
}
#endif  

#endif
