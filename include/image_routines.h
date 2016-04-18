/* Structure to hold an image */
/* FITS keywords indicated in brackets */
/* char pointers are to dynamically allocated strings */

#ifndef IMAGE_ROUTINES
#define IMAGE_ROUTINES

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fitsio.h"
#include "mathconst.h"
#include "newstring.h"
#include "ra_dec_string.h"

#ifdef __cplusplus
extern "C" {
#endif  

#define MAXDIM 7
#define IMTYPE 54321

enum PROJ {UNDEF, SIN, TAN, ARC, CAR};

struct image {
  int type;             /* Magic value to test for valid image */
  int bitpix;           /* Pixel numerical format in file (BITPIX) */
  int naxis;            /* Number of dimensions (NAXIS) */
  long size[MAXDIM];    /* Number of pixels in each dimension (NAXISn) */
  int niter;            /* Number of CLEAN iterations (NITER) */
  char *filename;       /* File name (--) */
  char *object;         /* Object name (OBJECT) */
  char *telescope;      /* Telescope name (TELESCOP) */
  char *instrument;     /* Instrument name (INSTRUME)*/
  char *observer;       /* Observer name (OBSERVER) */
  char *date_obs;       /* Date of observation (DATE-OBS) */
  char *timesys;        /* Time system (TIMESYS) */
  char *date_map;       /* Date of image creation (DATE or DATE-MAP) */
  char *bunit;          /* Image pixel units (BUNIT) */
  char *map_type;       /* Type of map (--) */
  char *origin;         /* Creating program (ORIGIN) */
  char *stokes;         /* Stokes parameter of image */
  char *ctype[MAXDIM];  /* Axis label (CTYPEn) */
  double crval[MAXDIM]; /* Axis value at reference point (CRVALn) */
  double crpix[MAXDIM]; /* Reference point pixel location (CRPIXn) */
  double cdelt[MAXDIM]; /* Axis value increment per pixel (CDELTn) */
  double crota[MAXDIM]; /* "axis" rotation value (CROTAn) */
  double bmaj;          /* Restoring beam major axis (BMAJ) */
  double bmin;          /* Restoring beam minor axis (BMIN) */
  double bpa;           /* Beam position angle (BPA) */
  double obsra;         /* Pointing center RA [deg] (OBSRA) */
  double obsdec;        /* Pointing center Dec [deg] (OBSDEC) */
  double equinox;       /* Equinox of coordinates [year] (EQUINOX or EPOCH) */
  double freq;          /* Frequency, MHz (--) */
  double noise;         /* Map noise from Difmap's estimate (NOISE) */
  long  npixels;        /* Number of pixels (--) */
  float *pixels;        /* Pointer to pixel data */
  float datamin;        /* Minimum pixel value */
  float datamax;        /* Maximum pixel value */
  enum PROJ projection; /* Projection type (code) */
};

/*--------------------------------------------------------------------*/

struct image *do_read(char *file);
void do_write_fits(struct image *mp, char *file);
void do_header(struct image *mp, char *file);
void delete_map(struct image *mp);
struct image *new_map(void);
void get_img_header(fitsfile *fptr, struct image *mp, int *status);
void set_img_header(fitsfile *fptr, struct image *mp, int *status);
void get_fits_string(fitsfile *fptr, char *keyword, char **param);
void set_fits_string(fitsfile *fptr, char *keyword, char **param);
void get_fits_double(fitsfile *fptr, char *keyword, double *param);
void set_fits_double(fitsfile *fptr, char *keyword, double *param);
void convert_date(char **date);
void copy_empty_map(struct image *p, struct image *src);
void copy_changed_map(struct image *p, struct image *src, float *p1);
void resize_map(struct image *p, int n, int m, int cx, int cy);

#ifdef __cplusplus
}
#endif

#endif
