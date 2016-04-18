/*--------------------------------------------------------------------
 * «mockcbiRoutines»
 * Libreria para correr mockcbi como subrutina
 * Suponemos que la polarizacion es LL. En el caso de una distinta, 
 * las mediciones en samples[].rdata[] pueden ubicarse en otras celdas.
 *--------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mockcbiRoutines.h"

#define SQR(x)  ((x)*(x))   // Funcion para calcular cuadrados.
#define PROGNAME  "mockcbiRoutines"

int debug = 0;

/*--------------------------------------------------------------------
 * Subrutina que calcula la atenuacion debida al haz primario.
 *--------------------------------------------------------------------*/

double ***atenuacion(struct image *fg_image, struct uvf_header **header_obs, 
		     int n_archivos, struct pbeam beam)
{
  int i, j, ix, iy;
  long k, t;    // k recorre xi
  double dx, dy, x, y;
  double lobs, mobs, obsra, obsdec, raimage, decimage, ra_offset = 0;
  double freq, arc;      // variables para el haz primario
  double ***attenu;      // arreglo a retornar.
  int x0, y0;            // pixel en el centro del haz primario.
  float *pixs = (float *) malloc (fg_image->npixels * sizeof (float));

  printf("Calculando atenuacion.\n");
  t = time(0);
  attenu = (double ***) malloc(n_archivos * sizeof(double**));

  raimage = fg_image->crval[0]*RPDEG;
  decimage = fg_image->crval[1]*RPDEG;

  for (i = 0; i < n_archivos; i++)
    {
      /*
      attenu[i] = (double **) malloc((fg_image->size[0] 
				      * fg_image->size[1] + 1) 
				     * sizeof(double*));
				     */
      attenu[i] = (double **) malloc((fg_image->npixels + 1) 
				     * sizeof(double*));
      obsra = ra_offset + header_obs[i]->obsra*RPDEG;
      obsdec = header_obs[i]->obsdec*RPDEG;  
      if (debug)
      printf("aten: pointing center: %s, %s (%g arcmin from image center)\n",
	     ra_string(obsra), dec_string(obsdec),
	     slaDsep(raimage, decimage, obsra, obsdec)/RPARCM);
      //printf("obsra: %d, obsdec: %d\n", obsra, obsdec);
      direccos(obsra, obsdec, raimage, decimage, &lobs, &mobs);
      dx = fg_image->cdelt[0] * RPDEG;   // en radianes
      dy = fg_image->cdelt[1] * RPDEG;   // en radianes
      x0 = (fg_image->crpix[0]-1.0) + lobs/dx;
      y0 = (fg_image->crpix[1]-1.0) + mobs/dy;
      dx = fg_image->cdelt[0] * 60;   // en arcmin
      dy = fg_image->cdelt[1] * 60;   // en arcmin

      printf("crpix[0]: %d, crpix[1]: %d\n",
	     (int) fg_image->crpix[0] ,(int) fg_image->crpix[1]);
      printf("xo: %d, yo: %d\n", x0, y0);
      //  samples_mod = g_mockcbi(cmb_image, fg_image, header_obs, samples_obs);
      k = 1;
      for (iy = 0; iy < fg_image->size[1]; iy++) 
	{
	  y = (iy - y0) * dy;
	  for (ix = 0; ix < fg_image->size[0]; ix++) 
	    {
	      x = (ix - x0) * dx;
	      arc = RPARCM * sqrt(x * x + y * y); // radio en radianes 
	      attenu[i][k] = (double *) malloc(header_obs[i]->nif 
					    * sizeof(double));
	      
	      for(j = 0; j < header_obs[i]->nif; j++)
		{
		  freq = header_obs[i]->iffreq[j] * 1e-9; 
		  attenu[i][k][j] = primary_beam(arc, freq, &beam);
		}
	      k++;
	    }
	}
    }
  t = time(0) - t;
  printf("time atten: %d\n", (int)t);
    
  for (i = 0; i < fg_image->npixels; i++) {
    pixs[i] = fg_image->pixels[i];
    fg_image->pixels[i] = attenu[0][i + 1][0];
  }
  do_write_fits(fg_image, "!atenuacion.fits");
  for (i = 0; i < fg_image->npixels; i++) {
    fg_image->pixels[i] = pixs[i];
  }
  free (pixs);
  
  return attenu;
}

/*--------------------------------------------------------------------
 * Subrutina que corre mockcbi
 *--------------------------------------------------------------------*/

void mockcbi_sub(struct image *cmb_image, 
		 struct image *fg_image,
		 struct uvf_header *header, 
		 struct uvf_sample *samples,
		 struct uvf_sample *samples_salida,
		 struct pbeam beam)
{
  struct leakage_table  leakage = {0, 0, 0, NULL};
  double fg_spind = 0.0;

  Status status;
  if(!cmb_image)
  {
    printf("Error con cmb_image\n");
    exit(0);
  }
  copiar_uvf_samples(samples, samples_salida, header->nsamp);
  status = do_clear("clear", header, samples_salida);   /* Visibilidades = 0 */
  if (status == FAILURE)
    exit(1);

  status = do_add_cmb("add_cmb", 1.0, 0.0, header, samples_salida,
		      cmb_image, NULL, NULL, fg_image, fg_spind, &beam, &leakage);

  //status = do_add_cmb("add_cmb", 1.0, 0.0, header, samples_salida,
  //	      cmb_image, fg_image, 0.0, &beam, &leakage);

  if (status == FAILURE)
    exit(1);
}

/*--------------------------------------------------------------------
 * Subrutina que copia 2 arreglos de estructuras uvf_samples
 *--------------------------------------------------------------------*/

void copiar_uvf_samples(struct uvf_sample *sin, struct uvf_sample *sout, int n)
{
  int i, j;

  for(j = 0; j < n; j++)
  {
    sout[j].date = sin[j].date;
    sout[j].u = sin[j].u;
    sout[j].v = sin[j].v;
    sout[j].w = sin[j].w;
    sout[j].inttim = sin[j].inttim;
    for(i = 0; i < MAX_DAT; i++)
      sout[j].rdata[i] = sin[j].rdata[i];
    sout[j].baseline = sin[j].baseline;
  }
}

/*--------------------------------------------------------------------
 * Aplana la imagen a intensidad = valor.
 *--------------------------------------------------------------------*/

void aplanar(struct image *p, double valor)
{
  int i;
  for(i = 0; i < p->npixels; i++)
    {
      p->pixels[i] = valor;
    }
  p->datamax = valor;
  p->datamin = valor;
}

/*--------------------------------------------------------------------
 * Multiplica las visibilidades por factor.
 *--------------------------------------------------------------------*/

void multiplicar_vis(struct uvf_header *header,
		     struct uvf_sample *samples, 
		     double factor) {
  int i, k, nsamp, nchan;

  nsamp = header->nsamp;
  nchan = header->nchan;

  for (i = 0; i < nsamp; i++) {      /* Loop through samples */
    for (k = 0; k < nchan; k++) {    /* Loop through channels */
      samples[i].rdata[3*k]   *= factor;    /* real part */
      samples[i].rdata[3*k+1] *= factor;    /* imaginary part */
      samples[i].rdata[3*k+2] /= SQR(factor);    /* imaginary part */
    }
  }
}

/*--------------------------------------------------------------------
 * Zero all the visibilities (whether flagged or not), leaving weights
 * unchanged. All IFs and polarizations are processed (nchan slots).
 *--------------------------------------------------------------------*/

Status do_clear(char *command, struct uvf_header *header,
		     struct uvf_sample *vis)
{
  long i, nsamp;
  int k, nchan;
  
  if (header == NULL || vis == NULL) {
    printf("%s: You must read a datafile before zeroing visibilities!\n",
	   command);
    return FAILURE;
  }

  nsamp = header->nsamp;
  nchan = header->nchan;

  for (i = 0; i < nsamp; i++) {      /* Loop through samples */
    for (k = 0; k < nchan; k++) {    /* Loop through channels */
      vis[i].rdata[3*k]   = 0.0f;    /* Zero real part */
      vis[i].rdata[3*k+1] = 0.0f;    /* Zero imaginary part */
    }
  }
  //printf("%s: Set all visibilities to (0.0,0.0).\n", command);
  return SUCCESS;
}


/*--------------------------------------------------------------------
 * Read Fits header and data from file named "file" and return pointers
 * to allocated structures "header" and "samples". Return status.
 *--------------------------------------------------------------------*/

Status do_read_uvdata(char *file, struct uvf_header **header,
		   struct uvf_sample **samples)
{
  int status = 0, error, i;
  fitsfile *infptr;
  char err_text[FLEN_STATUS];

  /* Allocate memory for header */

  *header = malloc(sizeof(struct uvf_header));
  if (*header == NULL) {
    printf("Insufficient memory\n");
    return FAILURE;
  }
  
  /* Read the input file */

  fits_open_file(&infptr, file, READONLY, &status);
  if (status) {
    fits_get_errstatus(status, err_text);
    printf("Unable to read file \"%s\" (%s)\n",
	    file, err_text);
    return FAILURE;
  } else {
    printf("Reading file \"%s\"\n", file);
  }
  error = get_header(infptr, *header);
  if (error) {
    free(*header);
    *header = NULL;
    return FAILURE;
  }

  /* Allocate memory for data */
  /* Note: this is not very sensible; we allocate memory for the
     maximum number of IFs and polarizations! */

  printf("Reading %d uv points for %d polarizations and %d frequency channels\n",
	  (*header)->nsamp, (*header)->nstokes, (*header)->nif);
  *samples = malloc(((*header)->nsamp)*sizeof(struct uvf_sample));
  if (*samples == NULL) {
    printf("insufficient memory for %d uv points\n", (*header)->nsamp);
    return FAILURE;
  }

  /* Read the data */
  
  fits_movabs_hdu(infptr,  1, NULL, &status);
  for (i=0; i < (*header)->nsamp; i++) {
    error = get_sample(infptr, *header, i+1, &(*samples)[i]);
    if (error) {
    free(*header);
    *header = NULL;
    free(*samples);
    *samples = NULL;
    return FAILURE;
    }
  }

  /* Close the input file */

  fits_close_file(infptr, &status);
  if (status) {
    fits_get_errstatus(status, err_text);
    printf("\nError reading file \"%s\" (%s)\n",
	    file, err_text);
    return FAILURE;
  }
  return SUCCESS;
}

/*--------------------------------------------------------------------
 * Create a new FITS file using the supplied header and data, copying
 * tables from "infile".
 *--------------------------------------------------------------------*/

Status do_write(char *command, char *file, struct uvf_header *header,
		    struct uvf_sample *samples, char *infile)
{
  int status = 0, i, hdunum;
  long gcount;
  fitsfile *infptr, *outfptr;
  char err_text[FLEN_STATUS];
  double refdat;
  char *filename = file;

  if (header == NULL || samples == NULL || infile == NULL) {
    printf("%s: No data to write!\n", command);
    return FAILURE;
  }

  /* Make a copy of the file name minus the leading ! if any */
  if (*filename == '!') filename++; 

  /* Create the output file */

  fits_create_file(&outfptr, file, &status);
  if (status) {
    fits_get_errstatus(status, err_text);
    printf("%s: Unable to create file \"%s\" (%s)\n",
	    command, filename, err_text);
    return FAILURE;
  } else {
    printf("%s: Creating file \"%s\"\n", command, filename);
  }

  /* Write the header */

  gcount = header->nsamp;
  refdat = (int) header->start_jd;
  status = put_header(outfptr, header, gcount, refdat, PROGNAME);

  /* Write the data */

  for (i=0; i < gcount; i++) {
    if (status) break;
    status = put_sample(outfptr, header, i+1, refdat, &samples[i]);
  }
  if (status) {
    return FAILURE;
  }

  /* Open the input file so that we can copy from it */

  fits_open_file(&infptr, infile, READONLY, &status);
  if (status) {
    fits_get_errstatus(status, err_text);
    printf("%s: Unable to read file \"%s\" (%s)\n",
	   command, infile, err_text);
    return FAILURE;
  }

  /* Copy history */

  copy_history(infptr, outfptr, infile);

  /* Copy extensions */

  fits_get_num_hdus(infptr, &hdunum, &status);
  printf("%s: Copying %d extensions from \"%s\" to \"%s\"\n",
	 command, hdunum-1, infile, filename);
  for (i=1; i < hdunum; i++) {
    fits_movabs_hdu(infptr, i+1, NULL, &status);
    fits_copy_hdu(infptr, outfptr, 0, &status);
  }
  
  /* Close input file */

  status = 0;
  fits_close_file(infptr, &status);
  fits_report_error(stdout, status);
  
  /* Close output file */

  status = 0;
  fits_close_file(outfptr, &status);
  fits_report_error(stdout, status);
  
  return SUCCESS;
}

/*--------------------------------------------------------------------
 * Add sky data, generated by Fourier transforming the supplied CMB
 * and foreground maps after multiplying by the primary beam,
 * to the visibility data.
 *
 * Arguments:
 *    command         command name (for diagnostic messages)
 *    scalefactor     scale factor (multiplies both CMB and foreground)
 *    ra_offset       RA of pointing center minus RA of UV header
 *                    (e.g., for the TRAIL of LEAD-TRAIL observations)
 *    header          header parameters for UV file (error if null)
 *    vis             visibility data
 *    cmb_image       CMB I image (error if null)
 *    q_image         CMB Q image (may be null)
 *    u_image         CMB U image (may be null)
 *    fg_image        foreground image (may be null)
 *    fg_spind        foreground spectral index (alpha, I ~ nu^alpha)
 *    beam            primary beam
 *    leakage_table   antenna leakage parameters
 *--------------------------------------------------------------------*/

Status do_add_cmb(char *command, double scalefactor,
			 double ra_offset, struct uvf_header *header,
			 struct uvf_sample *vis,
			 struct image *cmb_image,
			 struct image *q_image,
			 struct image *u_image,
			 struct image *fg_image, double fg_spind,
			 struct pbeam *beam,
			 struct leakage_table *q)
{
  int iff, nif, stokes, nstokes;
  long i, j, nsamp, nx, ny, pixels;
  int bl, i1, i2;
  double u, v, freq;
  float *image, *qmage, *umage;
  unsigned long nn[3];
  int ndim=2;
  complex cvis, ivis, qvis, uvis, vvis;
  float ivisr, qvisr, uvisr, vvisr, ivisi, qvisi, uvisi, vvisi;
  complex d, d1, d2, dterm, phfac;
  double ra, dec, obsra, obsdec, raimage, decimage;
  double dx, dy, du, dv, t, maxvis;
  double lobs, mobs, xobs, yobs, lphs, mphs, xphs, yphs;
  double phi;
  Status status = SUCCESS;

  /* Check that the visibility data exist */

  if (header == NULL || vis == NULL) {
    printf("%s: you must read a visibility file first!\n", command);
    return FAILURE;
  }

  /* Check that the input image exists and doesn't have any funny
     rotation of the coordinate system */

  if (cmb_image == NULL) {
    printf("%s: no sky image has been supplied!\n", command);
    return FAILURE;
  }
  if (cmb_image->crota[0] != 0.0 || cmb_image->crota[1] != 0.0) {
    printf("%s: the supplied image is rotated relative to N-S\n",
	   command);
    return FAILURE;
  }

  /* If Q and/or U images have been supplied, check that they are
     compatible with the I image (same size and projection). */

  if (q_image && !check_fg(q_image, cmb_image)) {
    printf("%s: the Q image does not match the I image in"
	   " location, size, or projection\n", command);
    return FAILURE;
  }
  if (u_image && !check_fg(u_image, cmb_image)) {
    printf("%s: the U image does not match the I image in"
	   " location, size, or projection\n", command);
    return FAILURE;
  }

  /* If a foreground image has been supplied, check that it is
     compatible (same size and projection as the CMB image) */

  if (fg_image && !check_fg(fg_image, cmb_image)) {
    printf("%s: the foreground image does not match the CMB image in"
	   " location, size, or projection\n", command);
    return FAILURE;
  }

  /* Size of input image */

  nx = cmb_image->size[0];
  ny = cmb_image->size[1];
  pixels = nx*ny;

  /* Array dimensions for the FFT routine (note order!) */

  nn[0] = 0;  /* not used */
  nn[1] = ny; /* pixels in y */
  nn[2] = nx; /* pixels in x */

  /* Pixel size of image (in radians) */

  dx = RPDEG*cmb_image->cdelt[0];
  dy = RPDEG*cmb_image->cdelt[1];
  //printf ("dx = %g, dy = %g\n", dx, dy);
  //printf ("nx = %d, ny = %d\n", nx, ny);

  /* Cell size in uv plane (in wavelengths) */

  du = 1.0/(nx*dx);
  dv = 1.0/(ny*dy);

  /* Make temporary arrays to hold a complex copy of the sky maps */

  image = malloc((pixels*2+2)*sizeof(float));
  qmage = malloc((pixels*2+2)*sizeof(float));
  umage = malloc((pixels*2+2)*sizeof(float));
  if (image == 0 || qmage == 0 || umage == 0) {
    printf("%s: insufficient memory\n", command);
    free(image);
    free(qmage);
    free(umage);
    return FAILURE;
  }
  /* Add sentinels at the end of the array (just in case) */
  image[pixels*2] = -777.0f;
  image[pixels*2+1] = -777.0f;

  /* Get and report uv pointing center and phase center for this file */

  /* Phase center of visibility data, including offset */
  ra = ra_offset + header->ra*RPDEG;
  dec = header->dec*RPDEG;
  /* Pointing center (primary beam) of visibility data, including offset */
  obsra = ra_offset + header->obsra*RPDEG;
  obsdec = header->obsdec*RPDEG;
  /* Reference point in the input sky image */
  raimage = cmb_image->crval[0]*RPDEG;
  decimage = cmb_image->crval[1]*RPDEG;
  if (debug) {
    printf("%s: pointing center: %s, %s (%g arcmin from image center)\n",
	   command,  ra_string(obsra), dec_string(obsdec),
	   slaDsep(raimage, decimage, obsra, obsdec)/RPARCM);
    printf("%s: phase center: %s, %s\n", command,
	   ra_string(ra), dec_string(dec));
  }

  /* Find the direction cosines of the pointing center relative to
     the reference point of the image */

  direccos(obsra, obsdec, raimage, decimage, &lobs, &mobs);

  /* Find this point in the pixel coordinate system of the image (where
     pixel numbers are zero-based) */

  xobs = (cmb_image->crpix[0]-1.0) + lobs/dx;
  yobs = (cmb_image->crpix[1]-1.0) + mobs/dy;
  if (debug)
    printf("%s: pointing center is at pixel (%g,%g)\n",
	   command, xobs, yobs);

  /* Check that it lies in the image! */

  if (xobs < 0 || xobs > nx || yobs < 0 || yobs > ny) {
    printf("%s: pointing center (%g,%g) is outside the range of the image\n",
	   command, xobs, yobs);
    free(image);
    free(qmage);
    free(umage);
    return FAILURE;
  }

  /* Find the direction cosines of the phase center relative to
     the reference point of the image */

  direccos(ra, dec, raimage, decimage, &lphs, &mphs);

  /* Find this point in the pixel coordinate system of the image */

  xphs = (cmb_image->crpix[0]-1.0) + lphs/dx;
  yphs = (cmb_image->crpix[1]-1.0) + mphs/dy;
  if (debug)
    printf("%s: phase center is at pixel (%g,%g)\n",
	   command, xphs, yphs);

  /* Process one IF at a time */
  
  nif = header->nif;
  nstokes = header->nstokes;
  nsamp = header->nsamp;
#if DEBUG
  printf("P: %ld %d %d\n", nsamp, nif, nstokes);
#endif

  for (iff = 0; iff < nif; iff++) {
    freq = header->iffreq[iff];
    if (debug)
      printf("%s: channel %d (%g GHz):\n", command, iff, freq*1e-9);

    /* Make a copy of the CMB image multiplied by the primary beam for
     * this frequency centered at the pointing center, and scale from
     * K to Jy/pixel. */

    if (debug)
      printf("Doing apply_beam\n");
    apply_beam(image, cmb_image, fg_image, fg_spind, 
	       xobs, yobs, freq*1e-9, beam);
    if (u_image)
      apply_beam(umage, u_image, NULL, 0.0, 
		 xobs, yobs, freq*1e-9, beam);
    if (q_image)
      apply_beam(qmage, q_image, NULL, 0.0, 
		 xobs, yobs, freq*1e-9, beam);


#if DEBUG
    if (image[pixels*2] != -777.0f || image[pixels*2+1] != -777.0f) {
      printf("apply_beam crashed\n");
      return FAILURE;
    }
#endif

    /* Do FFT to put the image in the uv plane; note offset of pointer
     * to image for Numerical Recipes convention */

    if (debug)
      printf("Doing fourn\n"); 
    fourn(image-1, nn, ndim, 1);
    if (u_image) fourn(umage-1, nn, ndim, 1);
    if (q_image) fourn(qmage-1, nn, ndim, 1);

#if DEBUG
    if (image[pixels*2] != -777.0f || image[pixels*2+1] != -777.0f) {
      printf("fourn crashed\n");
      return FAILURE;
      }
#endif

    /* Phase-rotate the visibility data for a pointing at the pointing
     * center, instead of at the center of the sky image. */

    if (debug) 
      printf("Doing phase_rotate\n");
    phase_rotate(image, cmb_image->size[0], cmb_image->size[1], 
		 -xphs, -yphs);
    if (q_image) phase_rotate(qmage, cmb_image->size[0], cmb_image->size[1], 
		 -xphs, -yphs);
    if (u_image) phase_rotate(umage, cmb_image->size[0], cmb_image->size[1], 
		 -xphs, -yphs);
    
#if DEBUG
    if (image[pixels*2] != -777.0f || image[pixels*2+1] != -777.0f) {
      printf("phase_rotate crashed\n");
      return FAILURE;
    }
#endif

    /* Loop through samples */

    maxvis = 0.0;
    for (i = 0; i < nsamp; i++) {
      /* printf("sample %ld\n", i); */
      /* Get (u,v) in wavelengths */
      u = vis[i].u * freq;
      v = vis[i].v * freq;
      /* interpolate in FT of image to get I visibility at this point,
	 (ivisr, ivisi) */
      status = getvis(u/du, v/dv, image, cmb_image->size[0],
		      cmb_image->size[1], &ivisr, &ivisi);
      if (status == FAILURE) {
	printf (" (u, v) = (%g, %g), freq = %g\n", u * 1e-2, v * 1e-2, freq);
	break;
      }
      if (q_image) {
	status = getvis(u/du, v/dv, qmage, cmb_image->size[0],
			cmb_image->size[1], &qvisr, &qvisi);
      } else {
	qvisr = qvisi = 0.0;
      }
      if (u_image) {
	status = getvis(u/du, v/dv, umage, cmb_image->size[0],
			cmb_image->size[1], &uvisr, &uvisi);
      } else {
	uvisr = uvisi = 0.0;
      }
      /* Stokes V is always zero, but we include it in the expressions
       in case it ever is needed */
      vvisr = vvisi = 0.0;

      ivis = Complex(ivisr, ivisi);
      if ((t=Cabs(ivis)) > maxvis)
	maxvis = t;
      qvis = Complex(qvisr, qvisi);
      uvis = Complex(uvisr, uvisi);
      vvis = Complex(vvisr, vvisi);

      /* Get the antenna numbers */
      
      bl = vis[i].baseline;
      i1 = (bl/256) - 1;     /* first antenna, zero-based */
      i2 = bl%256 - 1;       /* second antenna, zero-based */

      /* Get the polarization leakage D terms */

      if (q->do_leakage) {
	d1 = Complex(q->leakage[i1][iff].real, q->leakage[i1][iff].imag);
	d2 = Complex(q->leakage[i2][iff].real, q->leakage[i2][iff].imag);
	d = Cadd(d1, Conjg(d2));

	/* Baseline angle: note that this refers the antenna
	   orientation to the baseline, and is different for
	   each baseline */
	
	phi = atan2(vis[i].v, vis[i].u);
	phfac = Complex(cos(2.0*phi),sin(2.0*phi));
      }

      /* Loop for stokes */
      for (stokes = 0; stokes < nstokes; stokes++) {
	int s = 3*(nstokes*iff + stokes);
	double weight = vis[i].rdata[s+2];
	if (weight == 0.0) continue;
	/* printf("Sample %ld if %d stokes %d=%s: %g\n",
	   i, iff, stokes, stokes_label(header->stokes[stokes]),
	   weight); */
	switch (header->stokes[stokes]) {
	case -4: /* LR = Q - iU */
	  cvis = Cadd(qvis, Cmul(Complex(0,-1), uvis));
	  /* Leakage */
	  if (q->do_leakage) {
	    dterm = RCmul(-1, Cmul(phfac, Cmul(d, ivis)));
	    cvis = Cadd(cvis, dterm);
	  }
	  break;
	case -3: /* RL = Q + iU */
	  cvis = Cadd(qvis, Cmul(Complex(0,1), uvis));
	  /* Leakage */
	  if (q->do_leakage) {
	    dterm = RCmul(-1, Cmul(Conjg(phfac), Cmul(d, ivis)));
	    cvis = Cadd(cvis, dterm);
	  }
	  break;
	case -2: /* LL = I - V */
	  cvis = Cadd(ivis, RCmul(-1, vvis));
	case -1: /* RR = I + V */
	  cvis = Cadd(ivis, vvis);
	  break;
	default:
	  /* Add nothing */
	  printf("Invalid stokes: %d\n", header->stokes[stokes]);
	  status = FAILURE;
	  break;
	}
	vis[i].rdata[s]   += cvis.r*scalefactor;
	vis[i].rdata[s+1] += cvis.i*scalefactor;
      }
    }
    //printf("%s: ch %d maximum I visib amp = %g Jy\n", command, iff,
    //maxvis*fabs(scalefactor));
  }

  // printf("%s: added CMB data to visibilities (factor %g)",
  //command, scalefactor);
  if (q->do_leakage)
    printf(", including leakage terms\n");
  //else
  //printf("\n");

  /* Free work array */

  free(image);
  free(qmage);
  free(umage);

  return status;
}

/*--------------------------------------------------------------------
 * Add the CMB and (optional) foreground images into a complex array,
 * (imaginary parts zero) multiplying by the primary beam.
 *
 * Output:
 *    image: the new complex pixel array, attenuated by primary beam
 * Input:
 *    sky             sky image to be copied
 *    fg              foreground image to be added (or null pointer)
 *    fg_spind        spectral index for foreground image
 *    x0, y0          pixel coordinates of the center of the primary beam
 *    freq            frequency in GHz
 *    beam            primary beam parameters
 *--------------------------------------------------------------------*/

void apply_beam(float *image, struct image *sky,
		struct image *fg, double fg_spind,
		double x0, double y0, double freq,
		struct pbeam *beam)
{
  long ix, iy, pix = 0;
  double atten, x, y, arc, dx, dy;
  double scale=0.0, fg_scale=0.0;
  double value, pixel, t, maxpix = 0.0;
  //struct image *A = do_read ("../gaussianas6.fits");
  //char nombre [255];

  /* Conversion factor for CMB, to convert K to Jy/pixel */

  dx = sky->cdelt[0]*60.0;  /* arcmin */
  dy = sky->cdelt[1]*60.0;  /* arcmin */
  pixel = fabs((dx*RPARCM)*(dy*RPARCM)); /* pixel area in sr */

  if (sky) {
    double xf = ((PLANCK_H * 1e9)/(BOLTZ_K * TCMB)) * freq;
                            /*  x = h nu / k Tcmb, freq in GHz */
    double g = SQR(xf)*exp(xf)/SQR(exp(xf)-1.0);
    scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED))*g*pixel*SQR(freq*1e9);
  }

  /* Conversion factor for foreground, to convert K to Jy/pixel.
     This is defined so that the conversion factor is the same
     as for CMB at a fixed reference frequency of 30 GHz.
  */

  if (fg) {
    /* Foreground scale factor: same as CMB at 30 GHz, but with
       power-law spectrum, S propto freq^(fg_spind) */
    double ref_freq = 30.0;
    double xf = ((PLANCK_H * 1e9)/(BOLTZ_K * TCMB)) * ref_freq;
    double g = SQR(xf)*exp(xf)/SQR(exp(xf)-1.0);
    fg_scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED))*g*pixel*SQR(ref_freq*1e9)
      * pow(freq/ref_freq, fg_spind);
  }
  if (debug)
    printf(" 1 K = %g (CMB), %g (fg %g) Jy/pixel\n", scale, fg_scale,
	   fg_spind);
  
  /* Run through the pixels of the map; find the distance of each from
     the pointing center, and multiply by the primary beam response
     for that radius. */

  for (iy = 0; iy < sky->size[1]; iy++) {
    y = (iy - y0)*dy;
    for (ix = 0; ix < sky->size[0]; ix++) {
      x = (ix - x0)*dx;
      arc = RPARCM*sqrt(x*x + y*y); /* radius in radians */
      /* Primary beam response */ 
      atten = primary_beam(arc, freq, beam);
      //A->pixels[pix] = atten;
      value = 0.0;
      if (sky)
	value += scale*sky->pixels[pix];
      if (fg)
	value += fg_scale*fg->pixels[pix];
      image[pix*2] = value*atten; /* Real part */
      image[pix*2 + 1] = 0.0f; /* Imaginary part */
      t = fabs(image[pix*2]);
      if (t > maxpix) maxpix = t;
      pix++;
    }
  }
  //sprintf (nombre, "!atenuacion_%g.fits", freq);
  //do_write_fits(A, nombre);
  if (debug)
    printf("Maximum pixel after pb attenuation = %g Jy/pixel\n", maxpix);
}

/*
 * Interpolate in the visibility array to find the visibility at (u,v);
 * return real and imaginary parts in (visr, visi).
 * The complex-visibility array has dimension 2*nx*ny.
 */

Status getvis(double u, double v, float *image, long nx, long ny,
	      float *visr, float *visi)
{
  long i1, i2, j1, j2;
  double du, dv;
  double v11, v12, v21, v22;

  /* The units of u and v are the cell size. Check that they lie in the
     array */
  if (fabs(u) > (nx/2)-1 || fabs(v) > (ny/2)-1) {
    printf("Error in getvis: u,v = %g,%g\n", u, v);
    return FAILURE;
  }

  /* Wrap negative frequencies to positive */
  if (u < -1e-6) {
    u = nx + u;
  }
  if (v < -1e-6) {
    v = ny + v;
  }
  /* Corners of bounding cell and offset within cell */
  i1 = u; i2 = (i1+1)%nx; du = u - i1;
  j1 = v; j2 = (j1+1)%ny; dv = v - j1;
  //i1 = ((int) (u)) % nx; i2 = (i1+1)%nx; du = u - i1; // GUILLE
  //j1 = ((int) (v)) % ny; j2 = (j1+1)%ny; dv = v - j1; // GUILLE
  //if (i1 < 0 || i1 > nx || j1 < 0 || j2 > ny) { 
  if (i1 < 0 || i1 >= nx || j1 < 0 || j1 >= ny) { // GUILLE 21/11/2005
    printf("Error in getvis: u,v = %g,%g, %ld,%ld, %ld,%ld\n", u, v, i1, i2, j1, j2);
    printf("                 u, v = (%g, %g), du = %g, dv = %g\n", 
	   u , v, du, dv); // GUILLE
    return FAILURE;
  }
  /* Bilinear interpolation: real part */
  v11 = image[2*i1 + 2*nx*j1]; /* [i1, j1] */
  v12 = image[2*i1 + 2*nx*j2]; /* [i1, j2] */
  v21 = image[2*i2 + 2*nx*j1]; /* [i2, j1] */
  v22 = image[2*i2 + 2*nx*j2]; /* [i2, j2] */
  *visr = (1-du)*(1-dv)*v11 + (1-du)*dv*v12 + du*(1-dv)*v21 + du*dv*v22;
  /* Bilinear interpolation: imaginary part */
  v11 = image[1 + 2*i1 + 2*nx*j1]; /* [i1, j1] */
  v12 = image[1 + 2*i1 + 2*nx*j2]; /* [i1, j2] */
  v21 = image[1 + 2*i2 + 2*nx*j1]; /* [i2, j1] */
  v22 = image[1 + 2*i2 + 2*nx*j2]; /* [i2, j2] */
  *visi = (1-du)*(1-dv)*v11 + (1-du)*dv*v12 + du*(1-dv)*v21 + du*dv*v22;

  return SUCCESS;
}

/*--------------------------------------------------------------------
 * Phase rotate the visibility data in "image" to refer phase to point
 * (x,y) instead of (0,0).
 * Multiply pixel (i,j) by exp(-2 pi i (x/ni + y/nj))
 *--------------------------------------------------------------------*/

void phase_rotate(float *image, long nx, long ny,
		  double x, double y)
{
  long i, j, pix;
  double u, v, phase;
  float re, im, c, s;
  double du = TWOPI*x/nx;
  double dv = TWOPI*y/ny;

  pix = 0;
  /* v (j) axis changes more slowly */
  for (j=0; j < ny; j++) {
    /* Spatial frequency for this j */
    if (j < ny/2)
      v = dv*j;
    else
      v = dv*(j-ny);
    /* u (i) axis changes faster */
    for (i=0; i < nx; i++) {
      /* Spatial frequency for this j */
      if (i < nx/2)
	u = du*i;
      else
	u = du*(i-nx);
      /* Phase */
      phase = u + v;
      c = cos(phase);
      s = sin(phase);
      /* Phase rotation */
      re = image[pix];
      im = image[pix+1];
      image[pix] = re*c - im*s;
      image[pix+1] = re*s + im*c;
      pix += 2; /* step by two (complex) */
    }
  }
}

/*--------------------------------------------------------------------
 * Convert (ra,dec) to direction cosines (l,m) relative to
 * phase-tracking center (ra0, dec0). All in radians.
 * Reference: Synthesis Imaging in Radio Astronomy II, p.388.
 *--------------------------------------------------------------------*/

void direccos(double ra, double dec, double ra0, double dec0,
	      double *l, double *m)
{
  double delta_ra = ra - ra0;
  double cosdec = cos(dec);
  
  *l = cosdec * sin(delta_ra);
  *m = sin(dec) * cos(dec0) - cosdec * sin(dec0) * cos(delta_ra);
}

/*--------------------------------------------------------------------
 * Check two images for compatibility. At present, we just check the
 * number and size of the pixels. Return non-zero if compatible.
 *--------------------------------------------------------------------*/

int check_fg(struct image *img1, struct image *img2)
{
  if (img1->size[0] == img2->size[0] &&
      img1->size[1] == img2->size[1] &&
      img1->cdelt[0] == img2->cdelt[0] &&
      img1->cdelt[1] == img2->cdelt[1])
    return 1;
  else
    return 0;
}

/*--------------------------------------------------------------------
 * scale input  visibilities
 *--------------------------------------------------------------------*/

void scale (struct uvf_header *header,             /* scale visibilities  */
	    struct uvf_sample *samples, float ascale) {
  int i, k, nsamp, nchan;

  nsamp = header->nsamp;
  nchan = header->nchan;


  //  printf ("scaling vis by %g \n",*scale);
  printf ("scaling vis by %g \n",ascale);
  
  for (i = 0; i < nsamp; i++) {      /* Loop through samples */
    for (k = 0; k < nchan; k++) {    /* Loop through channels */
      //if (!ch_list[k])
      //continue;
      if (samples[i].rdata[3*k+2]  != 0) {
	samples[i].rdata[3*k]  *=  ascale;   
 	samples[i].rdata[3*k+1] *=  ascale; 
      }    /* weights */
    }
  }
}

/* rotar visibilidades en PI y toma complejo conjugado de las visibilidades */
void rotar_vis(struct uvf_header *header,
	       struct uvf_sample *samples)
{
  int i, k, nsamp, nchan;

  nsamp = header->nsamp;
  nchan = header->nchan;
  
  for (i = 0; i < nsamp; i++) {      /* Loop through samples */
    /* rotate baseline in PI and take complex conjugate of the visibilities  */
    if (samples[i].u < 0) {
      samples[i].u *= -1;
      samples[i].v *= -1;
      for (k = 0; k < nchan; k++) 
	  {
	    /*
	    if (!ch_list[k])
	    {
		samples[i].rdata[3 * k + 1] = 0;
		continue;
	      }
	    */
	    samples[i].rdata[3 * k + 1] *= -1;
	  }
    }
  }
}

double Irms (struct uvf_header **header,
	     struct uvf_sample **samples, int n_archivos) {
  int arch, iff, samp;
  double suma = 0;
  for (arch = 0; arch < n_archivos; arch++) {
    for (iff = 0; iff < header[arch]->nif; iff++) {
      for (samp = 0; samp < header[arch]->nsamp; samp++) {
	  suma += samples[arch][samp].rdata[iff * 3 + 2];
      }
    }
  }
  return sqrt (1 / suma);
  
}

void residual_vis(struct uvf_header *header,
		  struct uvf_sample *samples,
		  struct uvf_sample *samples_mod,  
		  struct uvf_sample *samples_res) {
  int i, k, nsamp, nchan;
  
  nsamp = header->nsamp;
  nchan = header->nchan;
  

  for (i = 0; i < nsamp; i++) {      /* Loop through samples */
    samples_res[i].u   = samples[i].u  ;
    samples_res[i].v   = samples[i].v  ;
    samples_res[i].w   = samples[i].w  ;
    samples_res[i].inttim   = samples[i].inttim  ;
    samples_res[i].date   = samples[i].date  ;
    samples_res[i].baseline   = samples[i].baseline  ;

    for (k = 0; k < nchan; k++) {    /* Loop through channels */
//      if (!ch_list[k])
//	{
//	  samples_res[i].rdata[3*k]   = 0;
//	  samples_res[i].rdata[3*k+1] = 0;
//	  samples_res[i].rdata[3*k+2] = 0;
//	  continue;
//	}
      samples_res[i].rdata[3*k]   = samples[i].rdata[3*k] - samples_mod[i].rdata[3*k] ;    /* real part */
      samples_res[i].rdata[3*k+1] = samples[i].rdata[3*k+1] - samples_mod[i].rdata[3*k+1]   ;    /* imaginary part */
      samples_res[i].rdata[3*k+2] = samples[i].rdata[3*k+2];    /* weights */
      //      printf("%f ",samples_res[i].rdata[3*k]);
    }
  }
}


/* NO FUNCIONA */

void dirtyMap (struct image *im, struct uvf_header **header, struct uvf_sample **samples, 
	       int nVis, double ***atten) {
  int i, j, nx = im->size[0], ny = im->size[1], arch, samp, chan, cont = 0;
  double raimage, decimage, obsra, obsdec, lobs, mobs, dx, dy, x0, y0, x, y, kx, ra_offset = 0.0;
  double seno, coseno, suma, xf, g, fg_scale, pixel, ref_freq = 30.0e9;

  xf = (PLANCK_H * ref_freq) / (BOLTZ_K * TCMB);
  g = SQR(xf)*exp(xf)/SQR(exp(xf)-1.0);
  fg_scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * SQR(ref_freq);   /* K to Jy sr^-1 */
  pixel = fabs(im->cdelt[0] * RPDEG * im->cdelt[1] * RPDEG);    /* pixel solid angle in ster */
  fg_scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * pixel * SQR(ref_freq);    /* K to Jy pixel^-1 */


  /* infinitesimal steps */
  dx = im->cdelt[0] * RPDEG;   // radians
  dy = im->cdelt[1] * RPDEG;   // radians
  for (arch = 0; arch < nVis; arch++) {
    // phase center "absolute" coordinates in radians 
    obsra = ra_offset + header[arch]->obsra*RPDEG;
    obsdec = header[arch]->obsdec*RPDEG;  
    // direction cosines of the phase center in the image coordinate system
    direccos(obsra, obsdec, raimage, decimage, &lobs, &mobs);
    // Find the phase center in the pixel coordinate system of the image (where pixel numbers are zero-based)
    x0 = (im->crpix[0] - 1.0) + lobs / dx;
    y0 = (im->crpix[1] - 1.0) + mobs / dy;
    
    for (i = 0; i < nx; i++) {
      x = (i - x0) * dx;                                // radians
      for (j = 0; j < ny; j++) {
	y = (j - y0) * dy;                                    // radians 
	
	im->pixels[i + j* nx] = 0;
	cont = 0;
	suma = 0;
	for(samp = 0; samp < header[arch]->nsamp; samp++) {
	  for(chan = 0; chan < header[arch]->nif; chan++) {
	    if (samples[arch][samp].rdata[chan * 3 + 2] != 0) {
	      kx = (samples[arch][samp].u * x + samples[arch][samp].v * y) * header[arch]->iffreq[chan];
	      seno = sin(2 * PI * kx);
	      coseno = cos(2 * PI * kx);
	      
	      suma += (samples[arch][samp].rdata[chan * 3] * coseno +
		       samples[arch][samp].rdata[chan * 3 + 1] * seno); // atten[arch][i + nx * j + 1][chan];
	      cont++;
	    }
	    im->pixels[i + j* nx] += suma / fg_scale / cont;
	  }
	}
      }
    }
  }
  printf ("fg_scale = %g\n", fg_scale);
}

