#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fitsio.h"
#include "slalib.h"
#include "uvsubs.h"
#include "mathconst.h"
#include "ra_dec_string.h"

#define DEBUG  0

/* --------------------------------------------------------------------
 * Read header data from UVF file into supplied structure; return 0
 * for success, 1 for failure.
 * -------------------------------------------------------------------- */

int get_header(fitsfile *fptr, struct uvf_header *header)
{
  int status = 0;
  int error;
  int i, k, hdunum, hdutype, bitpix, naxis, group;
  long pcount, gcount, naxes[MAX_DIM];
  char extname[68];
  int extver;
  char object[FLEN_VALUE];
  char telescop[FLEN_VALUE];
  char instrume[FLEN_VALUE];
  char ctype[FLEN_VALUE];
  char datobs[FLEN_VALUE], bunit[FLEN_VALUE], radecsys[FLEN_VALUE];
  char keyctype[FLEN_KEYWORD], keycrval[FLEN_KEYWORD];
  char keycdelt[FLEN_KEYWORD], keycrpix[FLEN_KEYWORD];
  char observer[FLEN_KEYWORD];
  double crval, cdelt, crpix;
  double equinox, obsra, obsdec, ra, dec, freq, bw;
  double min_date, max_date;
  int stokes[MAX_STOKES], nstokes, delta_stokes;
  int nc, nif, nval;
  struct uvf_sample sample;

  /* Read the standard header keywords */

  fits_read_imghdr(fptr, MAX_DIM, NULL, &bitpix, &naxis,
		   naxes, &pcount, &gcount, NULL, &status);
  fprintf(stderr, "  bitpix = %d, gcount=%ld\n", bitpix, gcount);
  fits_report_error(stderr, status);

  /* Check for UVF */

  if (status || gcount <= 1  || naxes[0] != 0) {
    fprintf(stderr, "This is not UV-FITS file\n");
    return 1;
  }

  /* Look for the keywords that we know about */

  fits_read_key(fptr, TSTRING, "OBJECT",   object,   NULL, &status);
  if (status == 202) {status = 0; object[0] = '\0';}
  fits_read_key(fptr, TSTRING, "TELESCOP", telescop, NULL, &status);
  if (status == 202) {status = 0; telescop[0] = '\0';}
  fits_read_key(fptr, TSTRING, "INSTRUME", instrume, NULL, &status);
  if (status == 202) {status = 0; instrume[0] = '\0';}
  fits_read_key(fptr, TSTRING, "DATE-OBS", datobs,   NULL, &status);
  if (status == 202) {status = 0; datobs[0] = '\0';}
  fits_read_key(fptr, TSTRING, "BUNIT",    bunit,    NULL, &status);
  if (status == 202) {status = 0; bunit[0] = '\0';}
  fits_read_key(fptr, TSTRING, "RADECSYS", radecsys, NULL, &status);
  if (status == 202) {status = 0; radecsys[0] = '\0';}
  fits_read_key(fptr, TSTRING, "OBSERVER", observer, NULL, &status);
  if (status == 202) {status = 0; observer[0] = '\0';}
  fits_read_key(fptr, TDOUBLE, "EQUINOX",  &equinox, NULL, &status);
  if (status == 202) { /* EQUINOX not found; use EPOCH */
    status = 0;
    /* fprintf(stderr, "Warning: EQUINOX keyword missing\n"); */
    fits_read_key(fptr, TDOUBLE, "EPOCH",  &equinox, NULL, &status);
  }
  fits_read_key(fptr, TDOUBLE, "OBSRA",  &obsra, NULL, &status);
  if (status == 202) {status = 0; obsra = 0.0;}
  fits_read_key(fptr, TDOUBLE, "OBSDEC", &obsdec, NULL, &status);
  if (status == 202) {status = 0; obsdec = 0.0;}
  fprintf(stderr, "  %s,%s,%s,%s,%s,%s,%s,%.2f,%s,%s\n", 
	  object, telescop, instrume, observer, datobs, bunit, radecsys,
	  equinox, ra_string(obsra*RPDEG), dec_string(obsdec*RPDEG));

  /* Look at the axis descriptions */

  error = 0;
  nc = 0;
  nstokes = 0;
  nif = 0;
  ra = dec = freq = bw = 0.0;
  stokes[0] = 0;
  delta_stokes = 0;
  fprintf(stderr, "  Axes (%d): ", naxis-1);
  for (i=1; i<naxis; i++) {
    status = 0;
    fits_make_keyn("CTYPE", i+1, keyctype, &status);
    fits_make_keyn("CRVAL", i+1, keycrval, &status);
    fits_make_keyn("CDELT", i+1, keycdelt, &status);
    fits_make_keyn("CRPIX", i+1, keycrpix, &status);
    fits_read_key(fptr, TSTRING, keyctype, ctype, NULL, &status);
    fits_read_key(fptr, TDOUBLE, keycrval, &crval, NULL, &status);
    fits_read_key(fptr, TDOUBLE, keycdelt, &cdelt, NULL, &status);
    fits_read_key(fptr, TDOUBLE, keycrpix, &crpix, NULL, &status);
    if (status)
      fits_report_error(stderr, status);
    if (i>1) fprintf(stderr, ",");
    fprintf(stderr, "%s(%ld)", ctype, naxes[i]);
    if (strcmp(ctype, "RA") == 0) {
      ra = crval;
      error += (naxes[i] != 1);
    } else if (strcmp(ctype, "DEC") == 0) {
      dec = crval;
      error += (naxes[i] != 1);
    } else if (strcmp(ctype, "COMPLEX") == 0) {
      nc = naxes[i];
      error += (nc != 3);
      error += (crval != 1.0);
    } else if (strcmp(ctype, "STOKES") == 0) {
      nstokes = naxes[i];
      error += (nstokes < 1 || nstokes > MAX_STOKES);
      delta_stokes = cdelt;
      for (k=0; k < MAX_STOKES; k++) {
	stokes[k] = crval + (1+k-crpix)*cdelt;
      }
      error += (nc == 0); /* COMPLEX must be an earlier axis than STOKES */
    } else if (strcmp(ctype, "FREQ") == 0) {
      freq = crval;
      bw = cdelt;
      error += (naxes[i] != 1);
    } else if (strcmp(ctype, "IF") == 0) {
      nif = naxes[i];
      error += (nif < 1);
      error += (crval != 1.0);
      error += (nstokes == 0); /* STOKES must be an earlier axis than IF */
      error += (nc == 0); /* COMPLEX must be an earlier axis than IF */
    } else {
      error += 1;
    }
  }
  fprintf(stderr, "\n");
  if (error > 0) {
    fprintf(stderr, "This is not an acceptable UV-FITS file\n");
    return 1;
  }
  fprintf(stderr, "  RA=%s, Dec=%s, Stokes=",
	  ra_string(ra*RPDEG), dec_string(dec*RPDEG));
  for (k = 0; k < nstokes; k++)
    fprintf(stderr,"%s,", stokes_label(stokes[k]));
  fprintf(stderr, " Freq=%.6g, BW=%.6g\n", freq, bw);
  nval = nc*nif; /* Number of data values at each uv point */

  /* Look at the parameter descriptions */

  header->index_u = -1;
  header->index_v = -1;
  header->index_w = -1;
  header->index_baseline = -1;
  header->index_date1 = -1;
  header->index_date2 = -1;
  header->index_inttim = -1;
  fprintf(stderr, "  Parameters (%ld): ", pcount);
  for (i=0; i<pcount; i++) {
    status = 0;
    fits_make_keyn("PTYPE", i+1, keyctype, &status);
    fits_read_key(fptr, TSTRING, keyctype, ctype, NULL, &status);
    fits_make_keyn("PSCAL", i+1, keyctype, &status);
    fits_read_key(fptr, TDOUBLE, keyctype, &header->pscal[i], NULL, &status);
    fits_make_keyn("PZERO", i+1, keyctype, &status);
    fits_read_key(fptr, TDOUBLE, keyctype, &header->pzero[i], NULL, &status);
    if (status)
      fits_report_error(stderr, status);
    if (i>0) fprintf(stderr, ",");
    fprintf(stderr, "%s", ctype);
    if (strcmp("UU---SIN", ctype) == 0 ||
	strcmp("UU",       ctype) == 0)
      header->index_u = i;
    else if (strcmp("VV---SIN", ctype) == 0 ||
	     strcmp("VV",       ctype) == 0)
      header->index_v = i;
    else if (strcmp("WW---SIN", ctype) == 0 ||
	     strcmp("WW",       ctype) == 0)
      header->index_w = i;
    else if (strcmp("BASELINE", ctype) == 0)
      header->index_baseline = i;
    else if (strcmp("INTTIM", ctype) == 0)
      header->index_inttim = i; /* optional */
    else if (strcmp("DATE", ctype) == 0) {
      if (header->index_date1 == -1)
	header->index_date1 = i;
      else
	header->index_date2 = i;
    } else {
      fprintf(stderr, "\nIgnoring unknown parameter \"%s\"\n", ctype);
    }      
  }
  fprintf(stderr, "\n");

  if (header->index_u < 0 || header->index_v < 0 || header->index_w < 0 ||
      header->index_baseline < 0 || header->index_date1 < 0) {
    fprintf(stderr, "A required parameter is missing\n");
    return 1;
  }

  if (pcount > MAX_PAR) {
    fprintf(stderr, "Too many group parameters in UVF file\n");
    return 1;
  }
  if (DEBUG) {
    for (i=0; i < pcount; i++)
      fprintf(stderr,"%g %g\n", header->pscal[i], header->pzero[i]);
  }

  /* Save header */

  strcpy(header->object, object);
  strcpy(header->telescop, telescop);
  strcpy(header->instrume, instrume);
  strcpy(header->bunit, bunit);
  strcpy(header->radecsys, radecsys);
  strcpy(header->observer, observer);
  header->equinox = equinox;
  header->obsra = obsra;
  header->obsdec = obsdec;
  header->ra = ra;
  header->dec = dec;
  header->freq = freq;
  header->bw = bw;
  header->nstokes = nstokes;
  header->delta_stokes = delta_stokes;
  for (k = 0; k < nstokes; k++)
    header->stokes[k] = stokes[k];
  header->nif = nif;
  header->nchan = nif*nstokes;
  header->nsamp = gcount;
  header->pcount = pcount;
  header->n_antennas = 0;

  /* Read the data to find the time range */

  max_date = 0.0;
  min_date = 1e20;
  for (group=1; group <= gcount; group++) {
    get_sample(fptr, header, group, &sample); 
    if (sample.date < min_date) min_date = sample.date;
    if (sample.date > max_date) max_date = sample.date;
  }
  header->start_jd = min_date;
  header->end_jd = max_date;
  fprintf(stderr, "  Time range: MJD %.8f to %.8f\n", min_date, max_date);

  /* Find the extensions */

  fits_get_num_hdus(fptr, &hdunum, &status);
  for (i=1; i < hdunum; i++) {
    fits_movabs_hdu(fptr, i+1, &hdutype, &status);
    if (hdutype == BINARY_TBL || hdutype == ASCII_TBL) {
      fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status);
      fits_read_key(fptr, TINT, "EXTVER", &extver, NULL, &status);
      fprintf(stderr, "  Extension %d is %s version %d (%s table)\n", i,
	      extname, extver, hdutype == BINARY_TBL ? "binary" : "ascii");
      if (hdutype == BINARY_TBL && strcmp(extname, "AIPS AN") == 0)
	get_antable(fptr, i+1, header);
      if (hdutype == BINARY_TBL && strcmp(extname, "AIPS FQ") == 0)
	get_fqtable(fptr, i+1, header);
    } else {
      fprintf(stderr, "  Skipping unknown extension of type %d\n", hdutype);
    }
  }

  return status ? 1 : 0;
}

/* --------------------------------------------------------------------
 * Get information from an AIPS antenna table;
 * return 0 for success, 1 for failure.
 * -------------------------------------------------------------------- */

int get_antable(fitsfile *fptr, int hdu, struct uvf_header *header)
{
  int status = 0, i;
  int hdutype, col_name, col_xyz, type_name, type_xyz;
  long nrows, rpt_name, rpt_xyz, wid_name, wid_xyz;
  char *string[1];
  int anynull;

  fits_movabs_hdu(fptr, hdu, &hdutype, &status);
  fits_get_num_rows(fptr, &nrows,  &status);
  header->n_antennas = (int) nrows;
  fits_get_colnum(fptr, CASEINSEN, "ANNAME", &col_name, &status);
  fits_get_coltype(fptr, col_name, &type_name, &rpt_name, &wid_name,
		   &status);
  fits_get_colnum(fptr, CASEINSEN, "STABXYZ", &col_xyz, &status);
  fits_get_coltype(fptr, col_xyz, &type_xyz, &rpt_xyz, &wid_xyz,
		   &status);
  if (status || type_name != TSTRING ||
                type_xyz != TDOUBLE || rpt_xyz != 3)
    fprintf(stderr, "    Warning: Invalid AIPS antenna table: %d %ld %d %ld\n",
	    type_name, rpt_name, type_xyz, rpt_xyz);
  for (i=0; i < nrows; i++) {
    string[0] = header->antenna[i].name;
    fits_read_col(fptr, TSTRING, col_name, i+1, 1, 1,
		  "NONAME", string, &anynull, &status);
    fits_read_col(fptr, TDOUBLE, col_xyz, i+1, 1, 1,
		  NULL, &(header->antenna[i].x), NULL, &status);
    fits_read_col(fptr, TDOUBLE, col_xyz, i+1, 2, 1,
		  NULL, &(header->antenna[i].y), NULL, &status);
    fits_read_col(fptr, TDOUBLE, col_xyz, i+1, 3, 1,
		  NULL, &(header->antenna[i].z), NULL, &status);
    if (DEBUG)
      fprintf(stderr, "    %s %8.3f %8.3f %8.3f\n", header->antenna[i].name,
	      header->antenna[i].x, header->antenna[i].y, header->antenna[i].z);
  }
  fprintf(stderr, "    %d antennas ", (int) nrows);
  for (i=0; i < nrows; i++)
    fprintf(stderr, " %s", header->antenna[i].name);
  fprintf(stderr, "\n");
  return status ? 1 : 0;
}

/* --------------------------------------------------------------------
 * Get information from an AIPS frequency table;
 * return 0 for success, 1 for failure.
 * -------------------------------------------------------------------- */

int get_fqtable(fitsfile *fptr, int hdu, struct uvf_header *header)
{
  int status = 0, i;
  int hdutype, col_iff, type_iff;
  long nrows, rpt_iff, wid_iff;
  double f;

  fits_movabs_hdu(fptr, hdu, &hdutype, &status);
  fits_get_num_rows(fptr, &nrows,  &status);
  fits_get_colnum(fptr, CASEINSEN, "IF FREQ", &col_iff, &status);
  fits_get_coltype(fptr, col_iff, &type_iff, &rpt_iff, &wid_iff,
		   &status);
  if (status || type_iff != TDOUBLE || nrows != 1 || rpt_iff !=
      header->nif)
    fprintf(stderr, "    Warning: Invalid AIPS frequency table\n");
  for (i=0; i < rpt_iff; i++) {
    fits_read_col(fptr, TDOUBLE, col_iff, 1, i+1, 1,
		  NULL, &f, NULL, &status);
    header->iffreq[i] = header->freq + f;
  }
  fprintf(stderr, "    %d IFs ", header->nif);
  for (i=0; i < header->nif; i++)
    fprintf(stderr, " %.2f", header->iffreq[i]*1e-9);
  fprintf(stderr, "\n");
  return status ? 1 : 0;
}

/* --------------------------------------------------------------------
 * Read a single sample from a UVF file into supplied structure;
 * return 0 for success, 1 for failure.
 * -------------------------------------------------------------------- */

int get_sample(fitsfile *fptr, struct uvf_header *header,
	       int number, struct uvf_sample *sample)
{
  int status = 0, k;
  double pdata[MAX_PAR];

  /* Read the random parameters */
  fits_read_grppar_dbl(fptr, number, 1, header->pcount, pdata, &status);

  /* Read the visibility data */
  fits_read_img_flt(fptr, number, 1, 3*(header->nchan),
		    0.0, sample->rdata, NULL, &status);

  /* Interpret the random parameters */
  k = header->index_u;
  if (k >= 0)
    sample->u = pdata[k]*header->pscal[k] + header->pzero[k];
  k = header->index_v;
  if (k >= 0)
    sample->v = pdata[k]*header->pscal[k] + header->pzero[k];
  k = header->index_w;
  if (k >= 0)
    sample->w = pdata[k]*header->pscal[k] + header->pzero[k];
  k = header->index_baseline;
  if (k >= 0)
    sample->baseline = pdata[k]*header->pscal[k] + header->pzero[k];
  sample->date = -2400000.5;
  k = header->index_date1;
  if (k >= 0)
    sample->date += pdata[k]*header->pscal[k] + header->pzero[k];
  k = header->index_date2;
  if (k >= 0)
    sample->date += pdata[k]*header->pscal[k] + header->pzero[k];
  k = header->index_inttim;
  if (k >= 0)
    sample->inttim = pdata[k]*header->pscal[k] + header->pzero[k];
  return status ? 1 : 0;
}

/* --------------------------------------------------------------------
 * Write a single sample to a UVF file from supplied structure;
 * return 0 for success, 1 for failure.
 * INTTIM parameter is written if available.
 * -------------------------------------------------------------------- */

int put_sample(fitsfile *fptr, struct uvf_header *header,
	       int group, double refdat, struct uvf_sample *sample)
{
  int status = 0;
  int pcount = 6;
  double pdata[MAX_PAR];
  pdata[0] = sample->u;
  pdata[1] = sample->v;
  pdata[2] = sample->w;
  pdata[3] = sample->baseline;
  pdata[4] = (int) sample->date - refdat;
  pdata[5] = (sample->date - (int)sample->date)*86400.0;
  if (header->index_inttim >= 0) {
    pcount = 7;
    pdata[6] = sample->inttim;
  }
  if (DEBUG) {
    int n1, n2;
    n1 = (int) pdata[3] / 256;
    n2 = (int) pdata[3] % 256;
    fprintf(stdout, "%d: %8g %8g %8g %d-%d %8g %8g\n", group, pdata[0],
	    pdata[1], pdata[2], n1, n2, pdata[4], pdata[5]);
  }
  fits_write_grppar_dbl(fptr, group, 1, pcount, pdata, &status);
  fits_write_img_flt(fptr, group, 1, 3*(header->nchan),
		     sample->rdata, &status);
  return status ? 1 : 0;
}

/* --------------------------------------------------------------------
 * Compare headers; return 0 if the datasets are compatible and can be
 * combined, 1 otherwise.
 *
 * Input:
 *   h1, h2:  pointers to the two headers.
 *   f1, f2:  file numbers (for inclusion in diagnostic messages).
 * -------------------------------------------------------------------- */

int comp_headers(struct uvf_header *h1, struct uvf_header *h2,
		 int f1, int f2)
{
  int error = 0, k;

  /* Object name must match */
  if (strcmp(h1->object, h2->object) != 0) {
    fprintf(stderr, "ERROR: Object \"%s\" (file %d) does not match \"%s\" (file %d)\n",
	    h1->object, f1, h2->object, f2);
    error = 1;
  }
  /* Brightness unit must match */
  if (strcmp(h1->bunit, h2->bunit) != 0) {
    fprintf(stderr, "ERROR: Bunit \"%s\" (file %d) does not match \"%s\" (file %d)\n",
	    h1->bunit, f1, h2->bunit, f2);
    error = 1;
  }
  /* Polarizations must match */
  if (h1->nstokes != h2->nstokes) {
    fprintf(stderr, "ERROR: Number of polarizations %d (file %d) does not match %d (file %d)\n",
	    h1->nstokes, f1, h2->nstokes, f2);
    error = 1;
  } else {
    for (k=0; k < h1->nstokes; k++) {
      if (h1->stokes[k] != h2->stokes[k]) {
	fprintf(stderr, "ERROR: Polarization %d=%s (file %d) does not match %s (file %d)\n",
		k, stokes_label(h1->stokes[k]), f1,
		stokes_label(h2->stokes[k]), f2);
	error = 1;
      }
    }
  }
  /* RA and Dec must match within 1 arcsec */
  if (fabs(h1->ra - h2->ra) > 1.0/3600.0) {
    fprintf(stderr, "ERROR: RA %.6f (file %d) does not match %.6f (file %d)\n",
	    h1->ra, f1, h2->ra, f2);
    error = 1;
  }
  if (fabs(h1->dec - h2->dec) > 1.0/3600.0) {
    fprintf(stderr, "ERROR: DEC %.6f (file %d) does not match %.6f (file %d)\n",
	    h1->dec, f1, h2->dec, f2);
    error = 1;
  }
  /* Number of channels must match */
  if (h1->nif != h2->nif) {
    fprintf(stderr, "ERROR: Number of channels %d (file %d) does not match %d (file %d)\n",
	    h1->nif, f1, h2->nif, f2);
    error = 1;
  }
  /* Number of parameters must match */
  if (h1->pcount != h2->pcount) {
    fprintf(stderr, "ERROR: Number of \"random parameters\" %d (file %d) does not match %d (file %d)\n",
	    h1->pcount, f1, h2->pcount, f2);
    error = 1;
  }
  /* Antenna tables must match */
  if (h1->n_antennas != h2->n_antennas) {
    fprintf(stderr, "ERROR: Number of antennas %d (file %d) does not match %d (file %d)\n",
	    h1->n_antennas, f1, h2->n_antennas, f2);
    error = 1;
  } else {
    int bad = 0, i;
    for (i=0; i<h1->n_antennas; i++) {
      if (h1->antenna[i].x != h2->antenna[i].x) bad = 1;
      if (h1->antenna[i].y != h2->antenna[i].y) bad = 1;
      if (h1->antenna[i].z != h2->antenna[i].z) bad = 1;
      if (strcmp(h1->antenna[i].name, h2->antenna[i].name)) bad = 1;
    }
    if (bad) {
      error = 1;
      fprintf(stderr, "ERROR: Antenna tables in files %d and %d do not match\n", f1, f2);
    }
  }

  return error;
}

/* --------------------------------------------------------------------
 * Output UV FITS file header using information in supplied structure.
 * Return 0 if successful, != 0 if an error occurred.
 * -------------------------------------------------------------------- */

int put_header(fitsfile *fptr, struct uvf_header *header, int gcount,
	       double refdat, char *origin)
{
  int naxis = 7;
  long pcount = 6;
  long naxes[7] = {0, 3, 1, 1, 1, 1, 1};
  float bscale = 1.0, bzero = 0.0;
  int one = 1,  zero = 0;
  double dfac1, dfac2, dref1, dref2, frac;
  char str[80];
  int status = 0;
  int ier, year, month, day;
  char datobs[16];

  /* Create the DATOBS string */

  slaDjcl(refdat, &year, &month, &day, &frac, &ier);
  if (ier != 0)
    fputs("Invalid DATOBS for FITS file.\n", stderr);
  sprintf(datobs, "%4.4d-%2.2d-%2.2d", year, month, day);

  /* Check for presence of INTTIM */

  if (header->index_inttim >=0)
    pcount = 7;

  /* Write required keywords. */

  naxes[2] = header->nstokes;
  naxes[4] = header->nif;
  fits_write_grphdr(fptr, 1, FLOAT_IMG, naxis, naxes, pcount, gcount,
		    1, &status);    
  fits_write_key(fptr, TFLOAT, "BSCALE", &bscale, 0, &status);
  fits_write_key(fptr, TFLOAT, "BZERO",  &bzero, 0, &status);
  
  /* Additional keywords. */

  fits_write_key(fptr, TSTRING, "OBJECT",   header->object,
		 "Source name", &status);
  fits_write_key(fptr, TSTRING, "TELESCOP", header->telescop,
		 0, &status);
  fits_write_key(fptr, TSTRING, "INSTRUME", header->instrume,
		 0, &status);
  fits_write_key(fptr, TSTRING, "OBSERVER", header->observer,
		 0, &status);
  fits_write_key(fptr, TSTRING, "DATE-OBS", datobs,
		 0, &status);
  fits_write_key(fptr, TSTRING, "BUNIT",    header->bunit,   
		 0,  &status);
  fits_write_key(fptr, TSTRING, "RADECSYS", header->radecsys,
		 0, &status);
  fits_write_key(fptr, TDOUBLE,  "EQUINOX",  &header->equinox, 
		 "Equinox of RA/Dec",   &status);
  fits_write_key(fptr, TDOUBLE, "EPOCH",     &header->equinox,
		 "Alternate name for EQUINOX", &status);
  fits_write_key(fptr, TDOUBLE, "OBSRA",    &header->obsra,
		 "Antenna pointing RA", &status);
  fits_write_key(fptr, TDOUBLE, "OBSDEC",   &header->obsdec,
		 "Antenna pointing Dec", &status);

  /* FITS coordinate parameters */
  
  fits_write_key(fptr, TSTRING, "CTYPE2", "COMPLEX", "1=real, 2=imag, 3=weight", &status);
  fits_write_key(fptr, TINT,    "CRVAL2", &one,  0, &status);
  fits_write_key(fptr, TINT,    "CDELT2", &one,  0, &status);
  fits_write_key(fptr, TINT,    "CRPIX2", &one,  0, &status);

  fits_write_key(fptr, TSTRING, "CTYPE3", "STOKES", "Correlator: -1=RR, -2=LL, -3=RL, -4=LR",&status);
  fits_write_key(fptr, TINT,    "CRVAL3", &header->stokes,  0, &status);
  fits_write_key(fptr, TINT,    "CDELT3", &header->delta_stokes, 0, &status);
  fits_write_key(fptr, TINT,    "CRPIX3", &one,       0, &status);

  fits_write_key(fptr, TSTRING, "CTYPE4", "FREQ","Frequency, Hz", &status);
  fits_write_key(fptr, TDOUBLE, "CRVAL4", &header->freq, 0, &status);
  fits_write_key(fptr, TDOUBLE, "CDELT4", &header->bw,   0, &status);
  fits_write_key(fptr, TINT,    "CRPIX4", &one,  0, &status);

  fits_write_key(fptr, TSTRING, "CTYPE5", "IF", "IF number",&status);
  fits_write_key(fptr, TINT,    "CRVAL5", &one, 0, &status);
  fits_write_key(fptr, TINT,    "CDELT5", &one, 0, &status);
  fits_write_key(fptr, TINT,    "CRPIX5", &one, 0, &status);

  fits_write_key(fptr, TSTRING, "CTYPE6", "RA", "Right ascension, degrees", &status);
  fits_write_key(fptr, TDOUBLE, "CRVAL6", &header->ra, 0, &status);
  fits_write_key(fptr, TINT,    "CDELT6", &zero, 0, &status);
  fits_write_key(fptr, TINT,    "CRPIX6", &one,  0, &status);

  fits_write_key(fptr, TSTRING, "CTYPE7", "DEC", "Declination, degrees", &status);
  fits_write_key(fptr, TDOUBLE, "CRVAL7", &header->dec, 0, &status);
  fits_write_key(fptr, TINT,    "CDELT7", &zero, 0, &status);
  fits_write_key(fptr, TINT,    "CRPIX7", &one,  0, &status);

  /* FITS random parameters. */

  fits_write_key(fptr, TSTRING, "PTYPE1", "UU", "baseline u projection, seconds", &status);
  fits_write_key(fptr, TINT,    "PSCAL1", &one,  0, &status);
  fits_write_key(fptr, TINT,    "PZERO1", &zero, 0, &status);
  fits_write_key(fptr, TSTRING, "PTYPE2", "VV", "baseline v projection, seconds", &status);
  fits_write_key(fptr, TINT,    "PSCAL2", &one,  0, &status);
  fits_write_key(fptr, TINT,    "PZERO2", &zero, 0, &status);
  fits_write_key(fptr, TSTRING, "PTYPE3", "WW", "baseline w projection, seconds", &status);
  fits_write_key(fptr, TINT,    "PSCAL3", &one,  0, &status);
  fits_write_key(fptr, TINT,    "PZERO3", &zero, 0, &status);

  fits_write_key(fptr, TSTRING, "PTYPE4", "BASELINE", "256*ANT1 + ANT2", &status);
  fits_write_key(fptr, TINT,    "PSCAL4", &one,  0, &status);
  fits_write_key(fptr, TINT,    "PZERO4", &zero, 0, &status);

  dfac1 = 1.0;
  dref1 = refdat + 2400000.5; /* convert MJD to JD */
  fits_write_key(fptr, TSTRING, "PTYPE5", "DATE", "UTC Julian Date part 1", &status);
  fits_write_key(fptr, TDOUBLE, "PSCAL5", &dfac1, "Days", &status);
  fits_write_key(fptr, TDOUBLE, "PZERO5", &dref1,  0, &status);
  
  dfac2 = 1.0/86400.0;
  dref2 = 0.0;
  fits_write_key(fptr, TSTRING, "PTYPE6", "DATE", "UTC Julian Date part 2", &status);
  fits_write_key(fptr, TDOUBLE, "PSCAL6", &dfac2, "Days/86400 (sec)", &status);
  fits_write_key(fptr, TFLOAT,  "PZERO6", &dref2, 0, &status);
  if (pcount >= 7) {
    fits_write_key(fptr, TSTRING, "PTYPE7", "INTTIM", "Integration time (sec)", &status);
    fits_write_key(fptr, TINT,   "PSCAL7", &one, 0, &status);
    fits_write_key(fptr, TINT,   "PZERO7", &zero, 0, &status);
  }

  sprintf(str, "%.12s (run by %.12s on %.12s)", origin, 
	  getenv("USER"), getenv("HOST"));
  fits_write_key(fptr, TSTRING, "ORIGIN", str," ", &status);

  /* Date/time of file creation */

  fits_write_date(fptr, &status);
  fits_report_error(stderr, status);
  return status;
}

/* --------------------------------------------------------------------
 * Copy history from one FITS file to another.
 * Return 0 if successful, != 0 if an error occurred.
 * -------------------------------------------------------------------- */

int copy_history(fitsfile *infptr, fitsfile *outfptr, char *infile)
{
  int status = 0;
  int i, nkeys;
  char card[81];

  /* Position files at main HDU */

  fits_movabs_hdu(outfptr, 1, NULL, &status);
  fits_movabs_hdu(infptr,  1, NULL, &status);

  /* Find the number of cards in the input header */

  fits_get_hdrspace(infptr, &nkeys, NULL, &status);

  /* Write an optional separator in the output history */

  if (infile != NULL) {
    sprintf(card, "History from file %.40s:", infile);
    fits_write_history(outfptr, "-------------------------------------------------------------------", &status);
    fits_write_history(outfptr, card, &status);
  }

  /* Read all the input cards and copy those that are HISTORY */

  for (i = 0; i < nkeys; i++) {
    fits_read_record(infptr, i+1, card, &status);
    if (strncmp(card, "HISTORY", 7) == 0)
      fits_write_record(outfptr, card, &status);
  }

  return status;
}

/* --------------------------------------------------------------------
 * Return character code for an integer Stokes polarization code.
 * Ref: NRAO Going AIPS, Chapter 6
 * -------------------------------------------------------------------- */

char *stokes_label(int code)
{
  static char *stokes[] = {"YX", "XY", "YY", "XX",  /* -8 to -5 */
		     "LR", "RL", "LL", "RR",  /* -4 to -1 */
		     "?",                     /* 0 = undefined */
                     "I", "Q", "U", "V"};     /* 1 to 4 */
  if (code < -8 || code > 4)
    code = 0;
  return stokes[code + 8];
}

