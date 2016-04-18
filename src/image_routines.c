/* Changes
 *
 * 2001-Apr-19: new.
 * 2001-Sep-25: moved ra_string and dec_string to separate file.
 * 2003-Aug-07: define code for Stokes PI (pseudo-I).
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fitsio.h"
#include "image_routines.h"
#include "mathconst.h"
#include "newstring.h"
#include "ra_dec_string.h"

/*--------------------------------------------------------------------*/
/* Read an image from the specified file "file" and store it in a
   dynamically allocated image structure; return a pointer to the
   structure, or NULL in the event of failure */
/*--------------------------------------------------------------------*/

struct image *do_read(char *file)
{
  struct image *mp;
  fitsfile *fptr;
  int status;
  int i;
  long size;
  char err_text[FLEN_STATUS];
  float nulval = FLT_MIN; /* value for blanked pixels */
  int anynul;
  float *p;
  char *dir;

  /* Open the file */

  status = 0;
  fits_open_file(&fptr, file, READONLY, &status);
  if (status) {
    fits_get_errstatus(status, err_text);
    fprintf(stdout, "Error: unable to read file \"%s\"\n(%s)\n",
	    file, err_text);
    return NULL;
  } else {
    fprintf(stdout, "Reading file \"%s\"\n", file);
  }

  /* Allocate memory for the image (note: even if a later operation
     fails, it is not necessary to free this memory; it will be freed
     on the next read command or on program exit) */

  mp = new_map();
  if (mp == NULL)
    return NULL;

  /* Complete the file name by prepending the working directory if
     appropriate */

  if (*file == '/' || strncmp(file, "http:/", 6) == 0 ||
      strncmp(file, "ftp:/", 5) == 0)
    dir = NULL;
  else
    dir = getenv("PWD"); /* only works with C-shell! */
  if (dir) {
    int l = strlen(dir) + 1 + strlen(file) + 1;
    mp->filename = malloc(l);
    sprintf(mp->filename, "%s/%s", dir, file);
  } else {
    mp->filename = newstring(file);
  }

  /* Read the header */

  status = 0;
  get_img_header(fptr, mp, &status);

  /* Read the data */

  size = (mp->npixels)*sizeof(float);
  if (size > 0) {

    /* Allocate memory for the pixels */

    if ((mp->pixels = malloc(size)) == NULL) {
      fprintf(stdout, "Error: not enough memory available to load image\n");
      mp->npixels = 0;
    }

    /* Load the pixel data */

    fits_read_img(fptr, TFLOAT, 1, mp->npixels, &nulval, mp->pixels,
		  &anynul, &status);
    if (anynul)
      fprintf(stdout, "Image contains blanked pixels\n");

    /* Find max and min */

    p = mp->pixels;
    mp->datamax = -FLT_MAX;
    mp->datamin = FLT_MAX;
    for (i = 1; i < mp->npixels; i++) {
      if (*p != nulval) {
	if (*p < mp->datamin) mp->datamin = *p;
	if (*p > mp->datamax) mp->datamax = *p;
      }
      p++;
    }
  }

  /* Close the file */

  fits_close_file(fptr, &status);
  if (status) {
    fits_get_errstatus(status, err_text);
    fprintf(stdout, "Error: while reading file \"%s\" (%s)\n",
	    file, err_text);
    return NULL;
  }

  return mp;
}

/*--------------------------------------------------------------------*
 * Writes an image from the specified dynamically allocated image 
 * structure mp.
 *--------------------------------------------------------------------*/

void do_write_fits(struct image *mp, char *file)
{
  fitsfile *fptr;
  int status;
  long size;
  char err_text[FLEN_STATUS];
  char *dir;

  /* Open the file */

  status = 0;
  fits_create_file(&fptr, file, &status);
  if (status) {
    fits_get_errstatus(status, err_text);
    fprintf(stdout, "Error: unable to read file \"%s\"\n(%s)\n",
	    file, err_text);
  } else {
    fprintf(stdout, "Reading file \"%s\" for writting\n", file);
  }

  /* Complete the file name by prepending the working directory if
     appropriate */

  if (*file == '/' || strncmp(file, "http:/", 6) == 0 ||
      strncmp(file, "ftp:/", 5) == 0)
    dir = NULL;
  else
    dir = getenv("PWD"); /* only works with C-shell! */
  if (dir) {
    int l = strlen(dir) + 1 + strlen(file) + 1;
    mp->filename = malloc(l);
    sprintf(mp->filename, "%s/%s", dir, file);
  } else {
    mp->filename = newstring(file);
  }


  /* Read the header */

  status = 0;
  set_img_header(fptr, mp, &status);
  if (status) {
    fits_get_errstatus(status, err_text);
    fprintf(stdout, "Error: while setting image header in \"%s\" (%s)\n",
	    file, err_text);
  }
  
  /* Read the data */
  
  size = (mp->npixels)*sizeof(float);
  if (size > 0) {

    /* Save the pixel data */

    fits_write_img(fptr, TFLOAT, 1, mp->npixels, mp->pixels, &status);
    if (status) {
      fits_get_errstatus(status, err_text);
      fprintf(stdout, "Error: while writting image in file \"%s\" (%s)\n",
	      file, err_text);
    }
  }

  /* Close the file */

  fits_close_file(fptr, &status);
  if (status) {
    fits_get_errstatus(status, err_text);
    fprintf(stdout, "Error: while reading file \"%s\" (%s)\n",
	    file, err_text);
  }
}

/*--------------------------------------------------------------------
 * List the header of the supplied structure "mp" in the file
 * named "file" (or stdout if file==NULL)
 *--------------------------------------------------------------------*/

void do_header(struct image *mp, char *file)
{
  FILE *fp = stdout;
  const char *unset  = "(unset)";
  const char *utc = "UTC";
  int i;

  /* Check for valid image */

  if (mp == NULL || mp->type != IMTYPE) {
    fprintf(stdout, "Error: object is not an image\n");
    return;
  }

  /* Open output (listing) file if required */

  if (file != NULL) {
    fp = fopen(file, "a");
    if (!fp) {
      fprintf(stdout, "Error: can't open file %s\n", file);
      return;
    }
  }

  /* List header */

  fprintf(fp, "\nParameters of image \"%s\"\n\n",
	  mp->filename ? mp->filename : unset);

  fprintf(fp, "  Number of pixels:     %ld\n", mp->npixels);
  fprintf(fp, "  Data format (bitpix): %d\n", mp->bitpix);
  fprintf(fp, "  Number of axes:       %d\n", mp->naxis);
  for(i = 0; i < mp->naxis; i++) {
    fprintf(fp, "  Axis%2d:               %-8s %5ld %.2f %.2f %.2f %.2f\n",
	    i+1, mp->ctype[i] ? mp->ctype[i] : unset,
	    mp->size[i], mp->crval[i], mp->crpix[i],
	    mp->cdelt[i], mp->crota[i]);
  }
  fprintf(fp, "  Object name:          %s\n",
	  mp->object ? mp->object : unset);
  fprintf(fp, "  Telescope:            %s\n",
	  mp->telescope ? mp->telescope : unset);
  fprintf(fp, "  Instrument:           %s\n",
	  mp->instrument ? mp->instrument : unset);
  fprintf(fp, "  Observer:             %s\n",
	  mp->observer ? mp->observer : unset);
  fprintf(fp, "  Observation date:     %s\n",
	  mp->date_obs ? mp->date_obs : unset);
  fprintf(fp, "  Time system:          %s\n",
	  mp->timesys ? mp->timesys : utc);
  fprintf(fp, "  Creation date:        %s\n",
	  mp->date_map ? mp->date_map : unset);
  fprintf(fp, "  Origin:               %s\n",
	  mp->origin ? mp->origin : unset);
  fprintf(fp, "  Map projection        ");
  switch (mp->projection) {
  case SIN: fprintf(fp, "SIN (orthographic)\n"); break;
  case TAN: fprintf(fp, "TAN (gnomonic)\n"); break;
  case CAR: fprintf(fp, "CAR\n"); break;
  case ARC: fprintf(fp, "ARC\n"); break;
  default:  fprintf(fp, "(unset)\n"); break;
  }
  if (mp->projection != UNDEF) {
    fprintf(fp, "  Coordinate equinox:   %c%.2f\n",
	    mp->equinox >= 1984.0 ? 'J' : 'B', mp->equinox);
    fprintf(fp, "  RA, Dec:              %s %s\n", 
	    ra_string(mp->crval[0]*RPDEG),
	    dec_string(mp->crval[1]*RPDEG));
    if (mp->obsra != 0.0 && mp->obsdec != 0.0) {
      fprintf(fp, "  Pointing RA, Dec:     %s %s\n",
	      ra_string(mp->obsra*RPDEG),
	      dec_string(mp->obsdec*RPDEG));
    }
  }
  fprintf(fp, "  Frequency (Hz):       %g\n", mp->freq);
  fprintf(fp, "  Stokes parameter:     %s\n", mp->stokes ? mp->stokes : unset);
  fprintf(fp, "  Beam:                 %.3f × %.3f, pa %.1f\n", 
	  mp->bmaj, mp->bmin, mp->bpa);
  fprintf(fp, "  Estimated noise:      %g\n", mp->noise);
  fprintf(fp, "  Minimum pixel value:  %g\n", mp->datamin);
  fprintf(fp, "  Maximum pixel value:  %g\n", mp->datamax);
  fprintf(fp, "  Pixel units:          %s\n",
	  mp->bunit ? mp->bunit : unset);

  /* Close output (listing) file */

  if (fp != stdout && fclose(fp) == EOF)
    fprintf(stdout, "Error: while writing file %s\n", file);
}

/*--------------------------------------------------------------------*/
/* Create a new empty image structure and return a pointer to it
   (or NULL in the event of failure) */
/*--------------------------------------------------------------------*/

struct image *new_map(void)
{
  struct image *p;
  int i;

  p = malloc(sizeof(struct image));
  if (p != NULL) {
    p->type = IMTYPE;
    p->bitpix = 0;
    p->naxis = 0;
    p->filename = NULL;
    p->object = NULL;
    p->telescope = NULL;
    p->instrument = NULL;
    p->observer = NULL;
    p->date_obs = NULL;
    p->timesys = NULL;
    p->date_map = NULL;
    p->bunit = NULL;
    p->map_type = NULL;
    p->origin = NULL;
    for (i = 0; i < MAXDIM; i++) {
      p->size[i] = 0;
      p->ctype[i] = NULL;
      p->crval[i] = 0.0;
      p->crpix[i] = 0.0;
      p->cdelt[i] = 1.0;
      p->crota[i] = 0.0;
    }
    p->npixels = 0;
    p->pixels = NULL;
    p->niter = 0;
    p->bmaj = 0.0;
    p->bmin = 0.0;
    p->bpa = 0.0;
    p->equinox  = 2000.0;
    p->freq = 0.0;
    p->stokes = NULL;
    p->datamax = 0.0;
    p->datamin = 0.0;
    p->projection = UNDEF;
  }
  if (p == NULL) {
    fputs("Error in routine \"new_map\": memory allocation failure\n", stdout);
  }
  return p;
}

/*--------------------------------------------------------------------*/
/* Delete the image structure pointed to by "mp" */
/*--------------------------------------------------------------------*/

void delete_map(struct image *mp)
{
  int i;

  /* Do nothing if the pointer is null */

  if (mp == NULL)
    return;

  /* Check that the pointer is valid */

  if (mp->type != IMTYPE) {
    fprintf(stdout, "Error in routine \"delete_map\": pointer does not point ot an image\n");
    exit(EXIT_FAILURE);
  }

  /* Free the components */

  free(mp->filename);
  free(mp->object);
  free(mp->telescope);
  free(mp->instrument);
  free(mp->observer);
  free(mp->date_obs);
  free(mp->date_map);
  free(mp->bunit);
  free(mp->map_type);
  free(mp->origin);

  //GUILLE 14/06/2006
  free(mp->stokes);
  free(mp->timesys);

  for (i = 0; i < MAXDIM; i++)
    free(mp->ctype[i]);
  free(mp->pixels);

  /* Free the image structure */

  free(mp);
}

/*--------------------------------------------------------------------
 * Read the header of the fits file (pointer "fptr") and fill
 * in the parameters in the image structure (pointer "mp").
 *--------------------------------------------------------------------*/

void get_img_header(fitsfile *fptr, struct image *mp, int *status)
{
  long pcount, gcount;
  int i;
  char err_text[FLEN_STATUS];
  char key[FLEN_KEYWORD];

  /* The following defines the symbolic names of the Stokes parameters
     coded as -8, -7, ..., 4, as defined by "Going AIPS", chapter 6,
     section 6.4.6, and in the draft paper "Representations of world
     coordinates in FITS" by E. W. Greisen & M. Calabreta (~2000).

     -9 is PI (Difmap).

     Other values are undefined. */

  static char *stokes[] = {"PI", "YX", "XY", "YY", "XX", "LR", "RL", "LL",
			   "RR", "none", "I", "Q", "U", "V"};

  if (*status != 0)
    return;
 
  /* Read the standard header keywords */

  fits_read_imghdr(fptr, MAXDIM, NULL, &(mp->bitpix), &(mp->naxis),
		   mp->size, &pcount, &gcount, NULL, status);
  if (*status) {
    fits_get_errstatus(*status, err_text);
    fprintf(stdout, "Error: %s\n", err_text);
    return;
  }
  if (mp->naxis < 1) {
    mp->npixels = 0;
    fprintf(stdout, "Error: file is not an image\n");
  } else {
    mp->npixels = mp->size[0];
    for (i=1; i < mp->naxis; i++)
      mp->npixels *= mp->size[i];
    if (mp->npixels < 1 || pcount != 0 || gcount != 1)
      fprintf(stdout, "Error: not an image (probably a UV FITS file)\n");
  }

  /* Look for the axis keywords */

  for (i = 0; i < mp->naxis; i++) {
    *status = 0;
    fits_make_keyn("CTYPE", i+1, key, status);
    get_fits_string(fptr, key, &mp->ctype[i]);
    fits_make_keyn("CRVAL", i+1, key, status);
    get_fits_double(fptr, key, &mp->crval[i]);
    fits_make_keyn("CDELT", i+1, key, status);
    get_fits_double(fptr, key, &mp->cdelt[i]);
    fits_make_keyn("CRPIX", i+1, key, status);
    get_fits_double(fptr, key, &mp->crpix[i]);
    fits_make_keyn("CROTA", i+1, key, status);
    get_fits_double(fptr, key, &mp->crota[i]);
    if (mp->ctype[i] && strcmp("FREQ", mp->ctype[i]) == 0) {
      mp->freq = mp->crval[i];
    }
    if (mp->ctype[i] && strcmp("STOKES", mp->ctype[i]) == 0) {
      int t = mp->crval[i];
      free(mp->stokes);
      if (t < -9 || t > 4) 
	mp->stokes = newstring("unknown");
      else
	mp->stokes = newstring(stokes[t+9]);
    }
  } 

  /* Find projection */

  if (mp->ctype[0] == NULL || mp->ctype[1] == NULL) {
    mp->projection = UNDEF;
  } else  if ((strcmp(mp->ctype[0], "RA---SIN") == 0) &&
	      (strcmp(mp->ctype[1], "DEC--SIN") == 0)){
    mp->projection = SIN;
  } else if ((strcmp(mp->ctype[0], "RA---TAN") == 0) &&
	     (strcmp(mp->ctype[1], "DEC--TAN") == 0)){
    mp->projection = TAN;
  } else if ((strcmp(mp->ctype[0], "RA---ARC") == 0) &&
	     (strcmp(mp->ctype[1], "DEC--ARC") == 0)){
    mp->projection = ARC;
  } else if ((strcmp(mp->ctype[0], "RA---CAR") == 0) &&
	     (strcmp(mp->ctype[1], "DEC--CAR") == 0)){
    mp->projection = CAR;
  } else {
    mp->projection = UNDEF;
  } 

  /* Look for the other keywords that we know about */

  get_fits_string(fptr, "OBJECT", &mp->object);
  get_fits_string(fptr, "TELESCOP", &mp->telescope);
  get_fits_string(fptr, "INSTRUME", &mp->instrument);
  get_fits_string(fptr, "OBSERVER", &mp->observer);
  get_fits_string(fptr, "ORIGIN", &mp->origin);
  get_fits_string(fptr, "DATE-OBS", &mp->date_obs);
  get_fits_string(fptr, "TIMESYS", &mp->timesys);
  get_fits_string(fptr, "DATE", &mp->date_map);
  if (mp->date_map == NULL)
      get_fits_string(fptr, "DATE-MAP", &mp->date_map);
  get_fits_string(fptr, "BUNIT", &mp->bunit);
  get_fits_double(fptr, "BMAJ", &mp->bmaj);
  get_fits_double(fptr, "BMIN", &mp->bmin);
  get_fits_double(fptr, "BPA", &mp->bpa);
  get_fits_double(fptr, "OBSRA", &mp->obsra);
  get_fits_double(fptr, "OBSDEC", &mp->obsdec);
  get_fits_double(fptr, "NOISE", &mp->noise);
  get_fits_double(fptr, "EQUINOX", &mp->equinox);
  if (mp->equinox == 0.0) 
    get_fits_double(fptr, "EPOCH", &mp->equinox);

  /* Check syntax of dates, and convert to 4-digit year
     if necessary */

  if (mp->date_obs != NULL)
    convert_date(&mp->date_obs);
  if (mp->date_map != NULL)
    convert_date(&mp->date_map);
}

/*--------------------------------------------------------------------
 * Writes the header of the fits file (pointer "fptr") from
 * the parameters in the image structure (pointer "mp").
 *--------------------------------------------------------------------*/

void set_img_header(fitsfile *fptr, struct image *mp, int *status)
{
  int i;
  char err_text[FLEN_STATUS];
  char key[FLEN_KEYWORD];

  if (*status != 0)
    return;
 
  /* Read the standard header keywords */

  fits_write_imghdr(fptr, mp->bitpix, mp->naxis,
		   mp->size, status);
  if (*status) {
    fits_get_errstatus(*status, err_text);
    fprintf(stdout, "Error: %s\n", err_text);
    return;
  }

  /* Look for the axis keywords */

  for (i = 0; i < mp->naxis; i++) {
    *status = 0;
    fits_make_keyn("CTYPE", i+1, key, status);
    set_fits_string(fptr, key, &mp->ctype[i]);
    fits_make_keyn("CRVAL", i+1, key, status);
    set_fits_double(fptr, key, &mp->crval[i]);
    fits_make_keyn("CDELT", i+1, key, status);
    set_fits_double(fptr, key, &mp->cdelt[i]);
    fits_make_keyn("CRPIX", i+1, key, status);
    set_fits_double(fptr, key, &mp->crpix[i]);
    fits_make_keyn("CROTA", i+1, key, status);
    set_fits_double(fptr, key, &mp->crota[i]);
    if (mp->ctype[i] && strcmp("FREQ", mp->ctype[i]) == 0) {
      mp->freq = mp->crval[i];
    }
    /*if (mp->ctype[i] && strcmp("STOKES", mp->ctype[i]) == 0) {
      int t = mp->crval[i];
      free(mp->stokes);
      if (t < -8 || t > 4)
	mp->stokes = newstring("unknown");
      else
	mp->stokes = newstring(stokes[t+8]);
	}*/
  } 

  /* Look for the other keywords that we know about */

  set_fits_string(fptr, "OBJECT", &mp->object);
  set_fits_string(fptr, "TELESCOP", &mp->telescope);
  set_fits_string(fptr, "INSTRUME", &mp->instrument);
  set_fits_string(fptr, "OBSERVER", &mp->observer);
  set_fits_string(fptr, "ORIGIN", &mp->origin);
  set_fits_string(fptr, "DATE-OBS", &mp->date_obs);
  set_fits_string(fptr, "TIMESYS", &mp->timesys);
  set_fits_string(fptr, "DATE", &mp->date_map);
  set_fits_string(fptr, "DATE-MAP", &mp->date_map);
  set_fits_string(fptr, "BUNIT", &mp->bunit);
  set_fits_double(fptr, "BMAJ", &mp->bmaj);
  set_fits_double(fptr, "BMIN", &mp->bmin);
  set_fits_double(fptr, "BPA", &mp->bpa);
  set_fits_double(fptr, "OBSRA", &mp->obsra);
  set_fits_double(fptr, "OBSDEC", &mp->obsdec);
  set_fits_double(fptr, "NOISE", &mp->noise);
  set_fits_double(fptr, "EQUINOX", &mp->equinox);
  set_fits_double(fptr, "EPOCH", &mp->equinox);

  /* Check syntax of dates, and convert to 4-digit year
     if necessary */

  /*if (mp->date_obs != NULL)
    convert_date(&mp->date_obs);
  if (mp->date_map != NULL)
  convert_date(&mp->date_map);*/
}

/*--------------------------------------------------------------------*/
/* Get value of a named string keyword; no error message if the
   keyword is absent */
/*--------------------------------------------------------------------*/

void get_fits_string(fitsfile *fptr, char *keyword, char **param)
{
  int status = 0;
  char value[FLEN_VALUE] = "";

  fits_read_key(fptr, TSTRING, keyword, value, NULL, &status);
  if (status == 0)
    *param = newstring(value);
  else
    *param = NULL;
}

/*--------------------------------------------------------------------*/
/* Set value of a named string keyword; no error message if the
   keyword is absent */
/*--------------------------------------------------------------------*/

void set_fits_string(fitsfile *fptr, char *keyword, char **param)
{
  int status = 0;

  fits_write_key(fptr, TSTRING, keyword, newstring(*param), NULL, &status);
  if (status != 0)
   fprintf(stdout, "Error: no se puede guardar %s al archivo.\n", keyword);
}

/*--------------------------------------------------------------------*/
/* Get value of a named numeric keyword; no error message if the
   keyword is absent */
/*--------------------------------------------------------------------*/

void get_fits_double(fitsfile *fptr, char *keyword, double *param)
{
  int status = 0;

  fits_read_key(fptr, TDOUBLE, keyword, param, NULL, &status);
  if (status != 0)
    *param = 0.0;
}

/*--------------------------------------------------------------------*/
/* Set value of a named numeric keyword; no error message if the
   keyword is absent */
/*--------------------------------------------------------------------*/

void set_fits_double(fitsfile *fptr, char *keyword, double *param)
{
  int status = 0;

  fits_write_key(fptr, TDOUBLE, keyword, param, NULL, &status);
  if (status != 0)
    fprintf(stdout, "Error: no se puede guardar %s al archivo.\n", keyword);
}

/*--------------------------------------------------------------------*/
/* Convert an old-style FITS format date "dd/mm/yy" to new format
 * "yyyy-mm-dd" if necessary; the argument is a pointer to a
 * dynamically allocated string which is reallocated if necessary. */
/*--------------------------------------------------------------------*/

void convert_date(char *date[])
{
  int day, month, year;
  char newdate[11];

  if ((*date)[2] == '/' && (*date)[5] == '/' && (*date)[8] == '\0') {
    day = atoi(*date);
    month = atoi(*date+3);
    year = atoi(*date+6) + 1900;
    sprintf(newdate, "%4.4d-%2.2d-%2.2d", year, month, day);
    free(*date);
    *date = newstring(newdate);
  }
}

/* --------------------------------------------------------------------
 * Copy header from one image to another; set all pixels to zero.
 * -------------------------------------------------------------------- */

void copy_empty_map(struct image *p, struct image *src)
{
  int i;

  p->type = IMTYPE;
  p->bitpix = -32;
  p->naxis = src->naxis;
  p->filename = NULL;
  p->object = newstring(src->object);
  p->telescope = newstring(src->telescope);
  p->instrument = newstring(src->instrument);
  p->observer = newstring(src->observer);
  p->date_obs = newstring(src->date_obs);
  p->timesys = newstring(src->timesys);
  p->date_map = NULL;
  p->bunit = newstring(src->bunit);
  p->map_type = newstring(src->map_type);
  p->origin = newstring(src->origin);
  for (i = 0; i < p->naxis; i++) {
    p->size[i] = src->size[i];
    p->ctype[i] = newstring(src->ctype[i]);
    p->crval[i] = src->crval[i];
    p->crpix[i] = src->crpix[i];
    p->cdelt[i] = src->cdelt[i];
    p->crota[i] = src->crota[i];
  }
  p->npixels = src->npixels;
  p->pixels = NULL;
  p->niter = 0;
  p->bmaj = 0.0;
  p->bmin = 0.0;
  p->bpa = 0.0;
  p->equinox  = src->equinox;
  p->freq = src->freq;
  p->noise = 0.0;
  p->datamax = 0.0;
  p->datamin = 0.0;
  p->projection = src->projection;
  p->obsra = 0.0;
  p->obsdec = 0.0;

  if (p->pixels == NULL && (p->pixels = malloc((p->npixels)*sizeof(float))) == NULL){
    fprintf(stdout, "Error: not enough memory for new image\n");
    p->npixels = 0;
  }

  for (i=0; i < p->npixels; i++) {
    p->pixels[i] = 0.0f;
  }
}

/* --------------------------------------------------------------------
 * Copy one image to another setting p as its pixels.
 * -------------------------------------------------------------------- */

void copy_changed_map(struct image *p, struct image *src, float *p1)
{
  int i;

  p->type = src->type;
  p->bitpix = src->bitpix;
  p->naxis = src->naxis;
  p->filename = NULL;
  p->object = newstring(src->object);
  p->telescope = newstring(src->telescope);
  p->instrument = newstring(src->instrument);
  p->observer = newstring(src->observer);
  p->date_obs = newstring(src->date_obs);
  p->timesys = newstring(src->timesys);
  p->date_map = NULL;
  p->bunit = newstring(src->bunit);
  p->map_type = newstring(src->map_type);
  p->origin = newstring(src->origin);
  for (i = 0; i < p->naxis; i++) {
    p->size[i] = src->size[i];
    p->ctype[i] = newstring(src->ctype[i]);
    p->crval[i] = src->crval[i];
    p->crpix[i] = src->crpix[i];
    p->cdelt[i] = src->cdelt[i];
    p->crota[i] = src->crota[i];
  }
  p->npixels = src->npixels;
  p->pixels = NULL;
  p->niter = 0;
  p->bmaj = src->bmaj;
  p->bmin = src->bmin;
  p->bpa = src->bpa;
  p->equinox  = src->equinox;
  p->freq = src->freq;
  p->noise = src->noise;
  p->datamax = p1[0];
  p->datamin = p1[0];
  p->projection = src->projection;
  p->obsra = src->obsra;
  p->obsdec = src->obsdec;

  if (p->pixels == NULL && (p->pixels = malloc((p->npixels)*sizeof(float))) == NULL){
    fprintf(stdout, "Error: not enough memory for new image\n");
    p->npixels = 0;
  }

  for (i=0; i < p->npixels; i++) {
    p->pixels[i] = p1[i];
    if (p1[i] < p->datamin)
      p->datamin = p1[i];
    if (p1[i] > p->datamax)
      p->datamax = p1[i];
  }
}

/* --------------------------------------------------------------------
 * Change de size of a 2 dimension image to n x m leaving the center 
 * at pixel (cx, cy).
 * -------------------------------------------------------------------- */

void resize_map(struct image *p, int n, int m, int cx, int cy)
{
  int x1, y1, x2, y2, k = 0;
  float *pixs;

  if (p->naxis < 2)
  {
    fprintf(stdout, "Error: Can't change map size, wrong number of dimensions: %d.\n",
	    p->naxis);
    return;
  }

  if ((pixs = malloc((n * m)*sizeof(float))) == NULL){
    fprintf(stdout, "Error: not enough memory for new image\n");
    return;
  }

  p->datamin = p->pixels[0];
  p->datamax = p->pixels[0];

  // Recorremos la imagen y la transformamos. 
  // (x, y) coordenadas de la nueva imagen.
  
  for (y2 = 0; y2 < m; y2++) 
    for (x2 = 0; x2 < n; x2++) 
    {
      x1 = x2 - (n / 2 - cx);
      y1 = y2 - (m / 2 - cy);
      if ((x1 < 0) || (x1 > p->size[0] - 1) ||
	  (y1 < 0) || (y1 > p->size[1] - 1))
	pixs[k] = 0;  // Out of image border.
      else
      {
	pixs[k] = p->pixels[y1 * p->size[0] + x1];
      }
      if (pixs[k] < p->datamin)
	p->datamin = pixs[k];
      if (pixs[k] > p->datamax)
	p->datamax = pixs[k];
      k++;
    }

  // Change the image size and value of CRPIX.

  p->size[0] = n;
  p->size[1] = m;
  p->npixels = n * m;
  p->crpix[0] = p->crpix[0] - (cx - n / 2);
  p->crpix[1] = p->crpix[1] - (cy - m / 2);

  free (p->pixels);
  p->pixels = pixs;
}
