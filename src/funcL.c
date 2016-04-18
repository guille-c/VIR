#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "funcL.h"

funcL * newFuncL (char * nombreImagen, char **nombresVis, int nVis, 
		  int n, int entropia, double cuantaSize,
		  double expanded_lx, double expanded_ly, double beamSize) {
  int i, j;
  double xf, g, bmaj, bmin;//, pixel;
  Status status;
  funcL *fL = malloc (sizeof (funcL));
  if (fL == NULL) {
    fprintf (stderr, "ERROR en newFuncL, fL = NULL.\n");
    exit (1);
  }
  //double ascale = 0.70e-3; // L1622_simu
  //double ascale = 0.0118788; // punt
  //double ascale = 1; // L1622

  /*
  L->malla = newMallaVoronoi();
  for (i = 0; i < n; i++) {
    insertarSitio (fL->malla, (double) random() / RAND_MAX,
		   (double) random() / RAND_MAX, 0.0);
  }
  */
  fL->malla = NULL;
  fL->nombreFits = malloc (50 * sizeof(char));
  if (fL->nombreFits == NULL) {
    fprintf (stderr, "ERROR en newFuncL, fL->nombre = NULL.\n");
    exit (1);
  }

  strcpy (fL->nombreFits, nombreImagen);
  fL->n_pols = n;
  fL->entropia = entropia;
  fL->ref_freq = 30.0e9; // 31.5e9;
  fL->expanded_lx = expanded_lx;
  fL->expanded_ly = expanded_ly;

  fL->chi2 = 0;
  fL->S = 0;
  fL->imagen = NULL;
  fL->cmb_image = do_read (nombreImagen);
  fL->fg_image = do_read (nombreImagen);
  fL->header_obs = NULL;
  fL->samples_obs = NULL;
  fL->samples_mod = NULL;
  //fL->Acutoff = cutoff;
  init_beam (&(fL->beam));

  if (beamSize == -2) {
    fL->beam.type = CBI2;  
    printf ("Usando beam CBI2.\n");
  }
  else if (beamSize <= 0) {
    fL->beam.type = CBI;
    printf ("Usando beam CBI.\n");
  }
  else {
    fL->beam.type = GAUSS;
    fL->beam.fwhm = beamSize * RPARCM;
    printf ("Usando beam GAUSS de %g arcmin.\n", beamSize);
  }
  // beam.cutoff *= 0.1;

  if (!fL->cmb_image) {
    fprintf (stderr, "ERROR al intentar abrir %s.", nombreImagen);
    exit (1);
  }

  if (!fL->cmb_image->bunit || !(strcmp("K", fL->cmb_image->bunit) == 0)) {
    fprintf (stderr, "ERROR: Unidades de %s  no estan en K. bunit = %s\n", 
	     nombreImagen, fL->cmb_image->bunit);
    //exit (1);
  }
  if (!fL->fg_image) {
    fprintf (stderr, "ERROR al intentar abrir %s.", nombreImagen);
    exit (1);
  }

  for (i = 0; i < fL->cmb_image->npixels; i++) {
    fL->cmb_image->pixels[i] = 0;
    fL->fg_image->pixels[i] = 0;
  }

  fL->n_archivos = nVis;
  fL->header_obs = (struct uvf_header**) malloc(fL->n_archivos * 
					    sizeof(struct uvf_header *));
  if (fL->header_obs == NULL) {
    fprintf (stderr, "ERROR en newFuncL, fL->header_obs = NULL.\n");
    exit (1);
  }

  fL->samples_obs = (struct uvf_sample **) malloc(fL->n_archivos * 
					     sizeof(struct uvf_sample *));
  if (fL->samples_obs == NULL) {
    fprintf (stderr, "ERROR en newFuncL, fL->samples_obs = NULL.\n");
    exit (1);
  }
  fL->samples_mod = (struct uvf_sample **) malloc(fL->n_archivos * 
					     sizeof(struct uvf_sample *));
  if (fL->samples_mod == NULL) {
    fprintf (stderr, "ERROR en newFuncL, fL->samples_mod = NULL.\n");
    exit (1);
  }
  //printf("ok2\n");
  fL->infile = malloc(fL->n_archivos * sizeof(char *));
  if (fL->infile == NULL) {
    fprintf (stderr, "ERROR en newFuncL, fL->infile = NULL.\n");
    exit (1);
  }
  for(i = 0; i < fL->n_archivos; i++) {
  // abrimos los archivos de visibilidades
    status = do_read_uvdata(nombresVis[i], &(fL->header_obs[i]), 
			    &(fL->samples_obs[i]));
    //multiplicar_vis(fL->header_obs[i], fL->samples_obs[i], norm_factor);
    if(status == SUCCESS) {
      fL->infile[i] = newstring(nombresVis[i]);
    }
    else {
      printf("Error con el archivo uvf\n");
      exit(1);
    }
  }
  //normalizarVisibilidades (fL->header_obs, fL->samples_obs, fL->n_archivos);
  for (i = 0; i < fL->n_archivos; i++) {
    fL->samples_mod[i] = (struct uvf_sample*) malloc((fL->header_obs[i]->nsamp) * 
						     sizeof(struct uvf_sample));
    if (fL->samples_mod[i] == NULL) {
      fprintf (stderr, "ERROR en newFuncL, fL->samples_mod[%d] = NULL.\n", i);
      exit (1);
    }
  }
  // ya tenemos abiertos los archivos de visibilidades

  xf = (PLANCK_H * fL->ref_freq) / (BOLTZ_K * TCMB);
  g = SQR(xf)*exp(xf)/SQR(exp(xf)-1.0);

  // ----------------RUIDO------------------------
  //pixel = fabs(fL->fg_image->cdelt[0] * RPDEG * fL->fg_image->cdelt[1] * RPDEG);
  /* pixel solid angle in ster */
  //L->fg_scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * pixel * SQR(fL->ref_freq);
  /* K to Jy pixel^-1 */
  fL->fg_scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * SQR(fL->ref_freq);
  /* K to Jy srad^-1 */

  bmaj = fL->fg_image->bmaj * RPDEG; // en radianes.
  bmin = fL->fg_image->bmin * RPDEG; // en radianes.
  printf ("bmaj = %g [deg] = %g [rad]\n", fL->fg_image->bmaj, bmaj);
  printf ("bmin = %g [deg] = %g [rad]\n", fL->fg_image->bmin, bmin);
  if (cuantaSize <= 0) {
    //fL->difmapNoise = 19.9451;  // [mJy / beam] IC443
    //fL->difmapNoise = 19.9451;  // [mJy / beam] punt
    //fL->difmapNoise = 11.0664;  // [mJy / beam] IC443
    //fL->difmapNoise = 9.29467;  // [mJy / beam] L1622 IRAS12
    //fL->difmapNoise = 16.8758;  // [mJy / beam] L1622 IRAS12
    //fL->difmapNoise = 7.41446;  // [mJy / beam] L1622_simu
    fL->difmapNoise = Irms (fL->header_obs, fL->samples_obs, fL->n_archivos);
    // Para considerar pixeles correlados
    double Nbeam = (PI * fL->fg_image->bmaj * fL->fg_image->bmin / (4 * log(2))/ 
		    fabs (fL->fg_image->cdelt[0] * fL->fg_image->cdelt[1]));
    printf ("Nbeam = %g\n", Nbeam);
    fL->difmapNoise *= sqrt (Nbeam);
  }
  else {
    fL->difmapNoise = cuantaSize;
  }
  printf ("quanta = %g [Jy / beam] \n", fL->difmapNoise);
  //fL->difmapNoise *= (4 * log(2) * 1e-3 / (PI * bmaj * bmin)); // [Jy / sr]
  fL->difmapNoise *= (4 * log(2) / (PI * bmaj * bmin)); // [Jy / sr]
  printf ("difmap noise 1 = %g [Jy / sr] \n", fL->difmapNoise);
  fL->difmapNoise /= fL->fg_scale; // en Kelvin.
  printf ("Imax = %g [Jy / beam] \n", 2.66 * (4 * log(2) / (PI * bmaj * bmin)) / fL->fg_scale);
  
  printf ("difmap noise 2 = %g K\n", fL->difmapNoise);
  printf ("fg_scale = %g\n", fL->fg_scale);

  /*
  for (arch = 0; arch < fL->n_archivos; arch++) {
    scale (fL->header_obs[arch], fL->samples_mod[arch], ascale);
    }*/
  /*
  printf ("Rotando visibilidades\n");
  for (arch = 0; arch < fL->n_archivos; arch++) {
    rotar_vis (fL->header_obs[arch], fL->samples_mod[arch]);
  }
  */
  /*
  {
    int arch, chan;
    for (arch = 0; arch < fL->n_archivos; arch++) {
      for (chan = 0; chan < fL->header_obs[arch]->nchan; chan++) {
	printf (" chan = %d, iffreq = %g\n", chan, 
		fL->header_obs[arch]->iffreq[chan]);
      }
    }
    exit (0);
  }
  */

  //do_write_fits (fL->fg_image, "!aux.fits");
  fL->atten = atenuacion(fL->fg_image, fL->header_obs, fL->n_archivos, fL->beam);

  fL->fourierI = (double **) malloc(fL->fg_image->size[0] * sizeof(double *));
  if (fL->fourierI == NULL) {
    fprintf (stderr, "ERROR en newFuncL, fL->fourierI = NULL.\n");
    exit (1);
  }
  fL->mask            = (int **) malloc(fL->fg_image->size[0] * sizeof(int *));
  if (fL->mask == NULL) {
    fprintf (stderr, "ERROR en newFuncL, fL->mask = NULL.\n");
    exit (1);
  }
  fL->PixIntegrados_I = (int **) malloc(fL->fg_image->size[0] * sizeof(int *));
  if (fL->PixIntegrados_I == NULL) {
    fprintf (stderr, "ERROR en newFuncL, fL->PixIntegrados_I = NULL.\n");
    exit (1);
  }
  fL->PixIntegrados_x = (int **) malloc(fL->fg_image->size[0] * sizeof(int *));
  if (fL->PixIntegrados_x == NULL) {
    fprintf (stderr, "ERROR en newFuncL, fL->PixIntegrados_x = NULL.\n");
    exit (1);
  }
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    fL->fourierI[i] = (double *) malloc (fL->fg_image->size[1]  * sizeof(double));
    if (fL->fourierI[i] == NULL) {
      fprintf (stderr, "ERROR en newFuncL, fL->fourier[%d] = NULL.\n", i);
      exit (1);
    }
    fL->mask[i]            = (int *) malloc (fL->fg_image->size[1]  * sizeof(int));
    if (fL->mask[i] == NULL) {
      fprintf (stderr, "ERROR en newFuncL, fL->mask[%d] = NULL.\n", i);
      exit (1);
    }
    fL->PixIntegrados_I[i] = (int *) malloc (fL->fg_image->size[1]  * sizeof(int));
    if (fL->PixIntegrados_I[i] == NULL) {
      fprintf (stderr, "ERROR en newFuncL, fL->PixIntegrados[%d] = NULL.\n", i);
      exit (1);
    }
    fL->PixIntegrados_x[i] = (int *) malloc (fL->fg_image->size[1]  * sizeof(int));
    if (fL->PixIntegrados_x[i] == NULL) {
      fprintf (stderr, "ERROR en newFuncL, fL->PixIntegrados_x[%d] = NULL.\n", i);
      exit (1);
    }
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      fL->PixIntegrados_I[i][j] = 0;
      fL->PixIntegrados_x[i][j] = 0;
    }
  }  

  printf ("Calculando sincos\n");
  calcularSin (fL);
  calculardxdy (fL);
  printf ("dx = %g, dy = %g\n", fL->dx[0] / RPARCM, fL->dy[0] / RPARCM);
  printf ("L listo\n");

  fL->nPixIntegrados_I = 0;
  printf ("nchan = %d, nif = %d\n\n", 
	  fL->header_obs[0]->nchan, fL->header_obs[0]->nif);
  return fL;
}

/* Calcula la funcion L.
 * En pars los iesimos parametros son la posicion x, los i+1 la posicion y
 * y los i+2 las intensidades de los n centros de Voronoi.
 */
double L (funcL *fL, double *pars, int n, int tipoMock) {
  int i, j, arch, samp, iff;
  double ret = 0, datamin = 1e300, datamax = -1e300, area = 0;
  NodoLista *nodo;
  //PoligonoVoronoi *pol;
  
  if (fL->malla != NULL) {
    eliminarMallaVoronoi(fL->malla);
    fL->malla = NULL;
  }
  
  fL->malla = newMallaVoronoi();
  for (i = 0; i < n; i++) {
    //printf ("ajustando (%g, %g, %g)\n", pars[3 * i], 
	    //pars[3 * i + 1], pars[3 * i + 2]* fL->difmapNoise);
    if (isnan (pars[3 * i]) || isnan (pars[3 * i + 1])) {
      printf ("ERROR en L, x[%d] = %g, y[%d] = %g\n", i, pars[3 * i], i, pars[3*i + 1]);
      //exit (1);
    }
    ajustarPoligonoCuadrado (&(pars[3 * i]), &(pars[3 * i + 1]));
    if (pars[3 * i + 2] <= 0) {
      pars[3 * i + 2] = MIN_PIX;
    }
    //printf ("insertando (%g, %g, %g -> %g)\n", pars[3 * i], 
    //pars[3 * i + 1], pars[3 * i + 2], pars[3 * i + 2] * fL->difmapNoise);
    insertarSitio (fL->malla, &(pars[3 * i]), &(pars[3 * i + 1]), 
		   pars[3 * i + 2] * fL->difmapNoise);

    if (datamin > pars[3 * i + 2] * fL->difmapNoise) {
      datamin = pars[3 * i + 2] * fL->difmapNoise;
    }
    if (datamax < pars[3 * i + 2] * fL->difmapNoise) {
      datamax = pars[3 * i + 2] * fL->difmapNoise;
    }
  }

  // Flujo
  /*
  for (i = 0, nodo = fL->malla->poligonos->primero; 
       i < n; i++, nodo = nodo->sgte) {
    pol = (PoligonoVoronoi *) nodo->info;
    pol->valor /= areaCuadrado (pol);
    area += areaCuadrado (pol);

    if (datamin > pol->valor) {
      datamin = pol->valor;
    }
    if (datamax < pol->valor) {
      datamax = pol->valor;
    }
  }
  */
  //printf ("--------AREA1 = %g\n", area);

  fL->fg_image->datamax = datamax;
  fL->fg_image->datamin = datamin;
  /*
  for (i = 0, nodo = fL->malla->poligonos->primero; 
       i < n; i++, nodo = nodo->sgte) {
    pol = (PoligonoVoronoi *) nodo->info;
    pol->valor = pars[n - 1 - i];
  }
  */
  if (fL->imagen != NULL) {
    for (i = 0; i < fL->fg_image->size[0]; i++) {
      free (fL->imagen[i]);
    }
    free (fL->imagen);
    fL->imagen = NULL;
  }
  fL->imagen = toImage (fL->malla, fL->fg_image->size[0], fL->fg_image->size[1], 
			fL->mask);
  
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      fL->fg_image->pixels[i + j * fL->fg_image->size[0]] = fL->imagen[i][j];
    }
  }

  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      if (fL->mask[i][j] == -1) {
	struct image *pix_image;

	fprintf(stderr, "ERROR en toImage: mask no esta completamente llena.\n");
	fprintf(stderr, "      (i, j) = (%d, %d).\n", i, j);
	imprimirMallaArchivo (fL->malla, "ERRORtoImage.dat");
	do_write_fits (fL->fg_image, "!ERRORtoImageImagen.fits");
	pix_image = do_read (fL->nombreFits);
	for (i = 0; i < pix_image->size[0]; i++) {
	  for (j = 0; j < pix_image->size[1]; j++) {
	    pix_image->pixels[i + j * pix_image->size[0]] =
	      fL->mask[i][j];
	  }
	}
	do_write_fits (pix_image, "!ERRORtoImage.fits");

	fL->imagen = toImageSinRaster (fL->malla, fL->fg_image->size[0], 
				       fL->fg_image->size[1], fL->mask);

	for (i = 0; i < pix_image->size[0]; i++) {
	  for (j = 0; j < pix_image->size[1]; j++) {
	    pix_image->pixels[i + j * pix_image->size[0]] =
	      fL->mask[i][j];
	  }
	}
	do_write_fits (pix_image, "!ERRORtoImageSinRaster.fits");

	for (i = 0; i < pix_image->size[0]; i++) {
	  for (j = 0; j < pix_image->size[1]; j++) {
	    pix_image->pixels[i + j * pix_image->size[0]] =
	      fL->imagen[i][j];
	  }
	}
	do_write_fits (pix_image, "!ERRORtoImageImagenSinRaster.fits");
	exit (1);
      }
    }
  }

  if (tipoMock == MOCK_EXACTO) {
    mockCBI(fL);
  }
  for (arch = 0; arch < fL->n_archivos; arch++) {
    if (tipoMock == MOCKCBI) {
      expanded_mockcbi(fL, arch);
//      mockcbi_sub (fL->cmb_image, fL->fg_image, 
//		   fL->header_obs[arch], fL->samples_obs[arch], 
//		   fL->samples_mod[arch], fL->beam);
    }
    for (iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
      for (samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	ret += ((SQR(fL->samples_mod[arch][samp].rdata[iff * 3] 
		     - fL->samples_obs[arch][samp].rdata[iff * 3])
		 + SQR(fL->samples_mod[arch][samp].rdata[iff * 3 + 1] 
		       - fL->samples_obs[arch][samp].rdata[iff * 3 + 1]))
		* fL->samples_mod[arch][samp].rdata[iff * 3 + 2]);
      }
    }
  }

  fL->chi2 = ret;
  if (fL->entropia) {
    fL->S = S (fL, pars, n);
    //printf ("L = chi2 + S = %g + %g = %g\n", fL->chi2 / 2, fL->S, fL->chi2 / 2 - fL->S);
    ret = fL->chi2 / 2 - fL->S;
  }
  else {
    //printf ("L = chi2 = %g\n", ret / 2);
    ret = fL->chi2 / 2;
  }
  //printf ("Chi2 = %g\n", ret);
  //free (im);
  /*  for (k = 0; k < fL->n_archivos; k++) {
    free (fL->samples_mod[k]);
    }*/
  //printf ("ret = %g\n", ret);

  return ret;
}

double S (funcL *fL, double *pars, int n) {
  int i;
  double S = 0, N = 0;
  
  for (i = 0; i < n; i++) {
    N += pars[3 * i + 2];
    S -= lgamma(pars[3 * i + 2] + 1);
  }
  S += lgamma (N + 1) - N * log (n);

  if (isinf(S)) {
    S = 0;
  }
  
  return S;
}

void dL (funcL *fL, double *pars, double *grad, int n, int aprox, int tipoMock){
  int i, j, t, arch;
  double datamin = 1e300, datamax = -1e300, area = 0;
  NodoLista *nodo;
  PoligonoVoronoi *pol;
  
  if (fL->malla != NULL){
    eliminarMallaVoronoi(fL->malla);
    fL->malla = NULL;
  }
  
  fL->malla = newMallaVoronoi();

  for (i = 0; i < n; i++) {
    ajustarPoligonoCuadrado (&(pars[3 * i]), &(pars[3 * i + 1]));
    
    if (pars[3 * i + 2] <= 0) {
      pars[3 * i + 2] = MIN_PIX;
    }

    insertarSitio (fL->malla, &(pars[3 * i]), &(pars[3 * i + 1]),
		   pars[3 * i + 2] * fL->difmapNoise);

    if (datamin > pars[3 * i + 2] * fL->difmapNoise) {
      datamin = pars[3 * i + 2] * fL->difmapNoise;
    }
    if (datamax < pars[3 * i + 2] * fL->difmapNoise) {
      datamax = pars[3 * i + 2] * fL->difmapNoise;
    }

    grad[3 * i] = 0;
    grad[3 * i + 1] = 0;
    grad[3 * i + 2] = 0;
  }

  // Flujo
  /*
  for (i = 0, nodo = fL->malla->poligonos->primero; 
       i < n; i++, nodo = nodo->sgte) {
    pol = (PoligonoVoronoi *) nodo->info;
    pol->valor /= areaCuadrado (pol);
    area += areaCuadrado (pol);

    if (datamin > pol->valor) {
      datamin = pol->valor;
    }
    if (datamax < pol->valor) {
      datamax = pol->valor;
    }
  }
  //printf ("--------AREA2 = %g\n", area);
  */

  fL->fg_image->datamax = datamax;
  fL->fg_image->datamin = datamin;
  
  if (fL->entropia) {
    dS(fL, pars, grad, n);
  }
  //printf ("Calculando ok2\n");

  fL->nPixIntegrados_I = 0;
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      fL->PixIntegrados_I[i][j] = 0;
      fL->PixIntegrados_x[i][j] = 0;
    }
  }

  fL->imagen = toImage(fL->malla, fL->fg_image->size[0], fL->fg_image->size[1], 
		       fL->mask);
  
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {

      //printf ("mask[%d][%d] = %d\n", i, j, fL->mask[i][j]);
      if (fL->mask[i][j] == -1) {
	struct image *pix_image;

	fprintf(stderr, "ERROR en toImage: mask no esta completamente llena\n");
	imprimirMallaArchivo (fL->malla, "ERRORtoImage.dat");
	pix_image = do_read (fL->nombreFits);
	for (i = 0; i < pix_image->size[0]; i++) {
	  for (j = 0; j < pix_image->size[1]; j++) {
	    pix_image->pixels[i + j * pix_image->size[0]] =
	      fL->mask[i][j];
	  }
	}
	do_write_fits (pix_image, "!ERRORtoImage.fits");	
	exit (1);
      }

      fL->fg_image->pixels[i + j * fL->fg_image->size[0]] = fL->imagen[i][j];
    }
  }

  if (tipoMock == MOCK_EXACTO) {
    mockCBI(fL);
  }
  else {
    for (arch = 0; arch < fL->n_archivos; arch++) {
      expanded_mockcbi(fL, arch);
//      mockcbi_sub(fL->cmb_image, fL->fg_image, 
//		  fL->header_obs[arch], fL->samples_obs[arch],
//		  fL->samples_mod[arch], fL->beam);
    }
  }

  guardarFits (fL->imagen, fL->mask, "dL", fL->nombreFits);


  printf ("Calculando dLdI\n");
  t = time (0);
  dLdI (fL, grad, n, aprox);
  printf ("tiempo dLdI = %d\n", (int) time (0) - t);
  printf ("Calculando dLdx\n");
  t = time(0);
  dLdx (fL, grad, n, aprox);
  printf ("tiempo dLdx = %d\n", (int) time (0) - t);
  printf ("Calculado dL\n");
  
  for (i = 0; i < n; i++) {
    printf ("grad[%d] = (%g, %g, %g)\n", i, grad[3*i], grad[3*i + 1], grad[3*i + 2]);
  }
  
  /*
  if (fL->malla != NULL){
    eliminarMallaVoronoi(fL->malla);
    fL->malla = NULL;
  }
  */
  // Liberamos la memoria de la imagen.
  if (fL->imagen != NULL) {
    for (i = 0; i < fL->fg_image->size[0]; i++) {
      free (fL->imagen[i]);
    }
    free (fL->imagen);
    fL->imagen = NULL;
  }
  else {
    printf ("NO ELIMINADA LA IMAGEN\n");
    exit(0);
  }
}

/* Calculo de menos el gradiente de la entropia. */

void dS (funcL *fL, double *pars, double *grad, int n){
  int i, j;
  double dS_cte = 0, N = 0;
  
  for (i = 0; i < n; i++) {
    N += pars[3 * i + 2];  
  }

  dS_cte = -log (n);
  for (i = 1; i <= N; i++) {
    dS_cte += 1.0 / i;
  }
  
  for (i = 0; i < n; i++) {
    grad[3 * i + 2] = -dS_cte;
    //printf ("N = %g, n = %d, N[%d] = %g\n", N, n, i, pars[3 * i + 2]);
    for (j = 1; j <= pars[3 * i + 2]; j++) {
      grad[3 * i + 2] += 1.0 / j;
    }
  }
}

/* Calculo del Gradiente con respecto a las posiciones de los
 * poligonos. 
 */

void dLdx (funcL *fL, double *grad, int n, int aprox){
  int i, arch, iff, samp;
  double dVdx_R, dVdx_I, dVdy_R, dVdy_I, Re, Im;
  PoligonoVoronoi *pol;
  NodoLista *nodo;

  // dLdx = 0;
/*   for (i = 0; i < n; i++) { */
/*     grad[3 * i]     = 0; */
/*     grad[3 * i + 1] = 0; */
/*   } */
/*   return; */

  for (arch = 0; arch < fL->n_archivos; arch++) {
    /*
    mockcbi_sub (fL->cmb_image, fL->fg_image, 
		 fL->header_obs[arch], fL->samples_obs[arch],
		 fL->samples_mod[arch], fL->beam);
    */
    /*for (i = 0, nodo = fL->malla->poligonos->primero; 
      i < 2; i++, nodo= nodo->sgte);*/
    for (i = 0, nodo = fL->malla->poligonos->primero; 
	 i < n; i++, nodo = nodo->sgte) {
      pol = (PoligonoVoronoi *) nodo->info;
      //printf ("pol = %d\n", pol->id);
      if (n - 1 - i != pol->id - 3) {
	fprintf (stderr, "ERROR: pol->id = %d != %d\n", 
		 pol->id - 3, n - 1 - i);
	exit (1);
      }
      //grad[3 * (n - 1 - i)]     = 0;
      //grad[3 * (n - 1 - i) + 1] = 0;
      Re = 0;
      Im = 0;
      for (iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
	for (samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	  dVdx_R = 0;
	  dVdx_I = 0;
	  dVdy_R = 0;
	  dVdy_I = 0;
	  dVdx (fL, arch, iff, samp, pol, &dVdx_R, &dVdx_I, &dVdy_R, &dVdy_I, aprox);
  	  Re = (fL->samples_mod[arch][samp].rdata[3 * iff] 
		- fL->samples_obs[arch][samp].rdata[3 * iff]) * dVdx_R;
	  Im = (fL->samples_mod[arch][samp].rdata[3 * iff + 1] 
		- fL->samples_obs[arch][samp].rdata[3 * iff + 1]) * dVdx_I;

	  grad[3 * (n - 1 - i)] += 
	    fL->samples_mod[arch][samp].rdata[3 * iff + 2] * (Re + Im);

	  if (isnan(grad[3 * (n - 1 - i)])) {
	    fprintf(stderr, "ERROR dLdx: grad[%d] = %g\n", 
		    3 * (n - 1 - i), grad[3 * (n - 1 - i)]);
	    fprintf(stderr, "          Re = %g, Im = %g\n", Re, Im);
	    fprintf(stderr, "          dVdxR = %g, dVdxI = %g\n", 
		    dVdx_R, dVdx_I);
	    fprintf(stderr, "          dVdyR = %g, dVdyI = %g\n", 
		    dVdy_R, dVdy_I);
	    fprintf(stderr, "          samp = %d, iff = %d\n", 
		    samp, iff);
	    exit (1);
	  }

	  Re = (fL->samples_mod[arch][samp].rdata[3 * iff] 
		- fL->samples_obs[arch][samp].rdata[3 * iff]) * dVdy_R;
	  Im = (fL->samples_mod[arch][samp].rdata[3 * iff + 1] 
		- fL->samples_obs[arch][samp].rdata[3 * iff + 1]) * dVdy_I;
	  grad[3 * (n - 1 - i) + 1] +=
	    fL->samples_mod[arch][samp].rdata[3 * iff + 2] * (Re + Im);

	  if (isnan(grad[3 * (n - 1 - i) + 1])) {
	    fprintf(stderr, "ERROR dLdy: grad[%d] = %g\n", 
		    3 * (n - 1 - i) + 1, grad[3 * (n - 1 - i) + 1]);
	    fprintf(stderr, "          Re = %g, Im = %g\n", Re, Im);
	    fprintf(stderr, "          dVdxR = %g, dVdxI = %g\n", 
		    dVdx_R, dVdx_I);
	    fprintf(stderr, "          dVdyR = %g, dVdyI = %g\n", 
		    dVdy_R, dVdy_I);
	    exit (1);
	  }
	}
      }
    }
  }
}

// Calculo del Gradiente con respecto a I. n es el numero de pols.

void dLdI (funcL *fL, double *grad, int n, int aprox) {
  int i, j, idPol;
  double *areasPix = (double *) malloc (n * sizeof(double)),
    *areas = (double *) malloc (n * sizeof(double)), area = 0;
  NodoLista *nodo;
  PoligonoVoronoi * pol;

  printf ("Calculando Fourier\n");
  if (aprox == APROX_DX) {
    calcularFourierIAprox (fL);
  }
  else if (aprox == APROX_VORONOI) {
    dLdI_Voronoi (fL, grad, n);
    return;
  }
  else {
    calcularFourierI (fL);
  }
  printf ("Calculado\n");

  // Inicializanmos el gradiente a 0.
  /*
  for (i = 0; i < n; i++) {
    grad[3 * i + 2] = 0;
  }
  */
//  // calculo areas==========================
//  for (i = 0; i < n; i++) {
//    areasPix[i] = 0;
//  }
//  for (i = 0; i < fL->fg_image->size[0]; i++) {
//    for (j = 0; j < fL->fg_image->size[1]; j++) {
//      idPol = fL->mask[i][j] - 3;
//      areasPix[idPol]++;
//    }
//  }
  nodo = fL->malla->poligonos->primero;
  for (i = 0; i < n; i++) {
    pol = (PoligonoVoronoi *) nodo->info;
    areas[pol->id - 3] =  areaCuadrado (pol);
    area += areas[pol->id - 3];
    nodo = nodo->sgte;
  }
  //========================================

  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      if (fL->mask[i][j] < 0) {
	fprintf (stderr, "ERROR en dLdI: fL->mask[%d][%d] = %d\n", 
		 i, j, fL->mask[i][j]);
	exit (1);
      }
      idPol = fL->mask[i][j] - 3;
      grad[3 * idPol + 2] += fL->fourierI[i][j] * fL->difmapNoise * fL->fg_scale;
      //grad[3 * idPol + 2] += fL->fourierI[i][j] * fL->difmapNoise * fL->fg_scale * areas[idPol] / areasPix[idPol]; Esto NO
      //grad[3 * idPol + 2] += fL->fourierI[i][j] * fL->difmapNoise * fL->fg_scale / areas[idPol] ;
    }
  }

  //Flujo
  /*
  for (i = 0, nodo = fL->malla->poligonos->primero; 
       i < n; i++, nodo = nodo->sgte) {
    pol = (PoligonoVoronoi *) nodo->info;
    printf ("calculando area para poligono {%g, %g, %g}....\n", pol->x, pol->y, pol->valor);
    grad[3 * pol->id + 2] /= areaCuadrado (pol);
    printf ("area = %g\n", areaCuadrado (pol));
  }
  */
  free (areas);
  free (areasPix);
  if (fL->nPixIntegrados_I != fL->fg_image->npixels) {
    struct image *pix_image, *mask;
    fprintf(stderr, "ERROR dL: Distinto numero de pixeles integrados\n");
    fprintf(stderr, "          nPixIntegrados = %d, npixels = %d\n",
	    fL->nPixIntegrados_I, (int) fL->fg_image->npixels);
    pix_image = do_read (fL->nombreFits);
    for (i = 0; i < pix_image->size[0]; i++) {
      for (j = 0; j < pix_image->size[1]; j++) {
	pix_image->pixels[i + j * pix_image->size[0]] =
	  fL->PixIntegrados_I[i][j];
      }
    }
    do_write_fits (pix_image, "!pixInt_ERROR.fits");
    
    mask = do_read (fL->nombreFits);
    for (i = 0; i < mask->size[0]; i++) {
      for (j = 0; j < mask->size[0]; j++) {
	mask->pixels[i + j * mask->size[0]] = fL->mask[i][j];
      }
    }
    do_write_fits(mask, "!mask_ERROR.fits");

    exit (1);
  }
}
  
void dLdI_Voronoi (funcL *fL, double *grad, int n) {
  int i, j, idPol, ix, iy, nx, ny, arch, iff, samp;
  double area = 0, x, y, u, v, kx, Re, Im;
  NodoLista *nodo;
  PoligonoVoronoi * pol;

  printf ("APROX VORONOI\n");
  // Inicializanmos el gradiente a 0.
  /*
  for (i = 0; i < n; i++) {
    grad[3 * i + 2] = 0;
  }
  */
  nx = fL->fg_image->size[0];
  ny = fL->fg_image->size[1];
  nodo = fL->malla->poligonos->primero;
  for (i = 0; i < n; i++) {
    pol = (PoligonoVoronoi *) nodo->info;
    ix = pol->x * (nx - 1);
    iy = pol->y * (ny - 1);
    
    area = areaCuadrado(pol) ;

    for (arch = 0; arch < fL->n_archivos; arch++) {
      y = (ix - fL->y0[arch]) * fL->dy[arch];
      x = (iy - fL->x0[arch]) * fL->dx[arch];
      for (iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
	for(samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	  u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
	  v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
	  kx = (u * x + v * y);
	  // Transf Fourier Directa
	  
	  Re = ((fL->samples_mod[arch][samp].rdata[3 * iff] 
		 - fL->samples_obs[arch][samp].rdata[3 * iff]) 
		* cos (2 * PI * kx));
	  Im = ((fL->samples_mod[arch][samp].rdata[3 * iff + 1] 
		 - fL->samples_obs[arch][samp].rdata[3 * iff + 1])
		* sin (2 * PI * kx));

	  //fL->fourierI[i][j] += fL->samples_mod[arch][samp].rdata[3 * iff + 2] 
	  //* (Re + Im);
	  grad [3 * (pol->id - 3) + 2] += fL->samples_mod[arch][samp].rdata[3 * iff + 2] 
	    * (Re + Im) * fL->atten[arch][ix + iy * fL->fg_image->size[0] + 1][iff];
	}
	//fL->fourierI[i][j] *= fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff];
      }
      grad [3 * (pol->id - 3) + 2] *= (area * (nx - 1) * (ny - 1) 
				       * fabs (fL->dx[arch] * fL->dy[arch])
				        * fL->difmapNoise * fL->fg_scale);
      printf ("grad[%d] = %g\n", 3 * (pol->id - 3) + 2, grad [3 * (pol->id - 3) + 2]);
    }
        
    nodo = nodo->sgte;
  }
}
  
void calcularSin (funcL *fL) {
  int arch, samp, iff;
  double dx, dy;//, x, y;
  double u, v, sinu, sinv;

  //dx = fL->fg_image->cdelt[0] * 60;   // en arcmin
  //dy = fL->fg_image->cdelt[1] * 60;   // en arcmin
  dx = fL->fg_image->cdelt[0] * RPDEG;   // en radianes
  dy = fL->fg_image->cdelt[1] * RPDEG;   // en radianes

  //printf ("dx = %g, dy = %g\n", dx, dy);

  fL->sin = (double ***) malloc (fL->n_archivos * sizeof (double **));
  if (fL->sin == NULL) {
    fprintf (stderr, "ERROR en calcularSin, fL->sin = NULL.\n");
    exit (1);
  }
  for (arch = 0; arch < fL->n_archivos; arch++) {
    fL->sin[arch] = (double **) malloc (fL->header_obs[arch]->nsamp 
					* sizeof(double*));
    if (fL->sin[arch] == NULL) {
      fprintf (stderr, "ERROR en calcularSin, fL->sin[%d] = NULL.\n", arch);
      exit (1);
    }
    
    for(samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
      fL->sin[arch][samp] = (double *) malloc(fL->header_obs[arch]->nif
					     * sizeof(double));
      if (fL->sin[arch][samp] == NULL) {
	fprintf (stderr, "ERROR en calcularSin, fL->sin[%d][%d] = NULL.\n", arch, samp);
	exit (1);
      }
      for (iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
	u = fL->samples_obs[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
	v = fL->samples_obs[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
	if (u != 0) {
	  sinu = sin (PI * u * fabs (dx)) / (PI * u);
	}
	else {
	  sinu = fabs (dx);
	  //sinu = dx / RPDEG;
	  printf ("sinu 0, samp = %d, iff = %d, u = %g, iffreq = %g\n", 
		  samp, iff, fL->samples_obs[arch][samp].u, 
		  fL->header_obs[arch]->iffreq[iff]);
	}
	if (v != 0) {
	  sinv = sin (PI * v * fabs (dy)) / (PI * v);
	}
	else {
	  sinv = fabs (dy);
	  //sinv = dy / RPDEG;
	  printf ("sinv 0, samp = %d, iff = %d, v = %g, iffreq = %g\n",
		  samp, iff, fL->samples_obs[arch][samp].u,
		  fL->header_obs[arch]->iffreq[iff]);
	}
	/*
	L->sin[arch][samp][iff] = (sin (PI * u * dx) * sin (PI * v * dy) 
				    / (PI * PI * u * v));
	*/
	fL->sin[arch][samp][iff] = sinu * sinv;
	if (isnan(fL->sin[arch][samp][iff])) {
	  fprintf(stderr, "ERROR calcularSin: sin[%d][%d][%d] = %g\n", 
		  arch, samp, iff, fL->sin[arch][samp][iff]);
	  exit (1);
	}
	//printf ("sinu = %g, sinv = %g, u = %g, v = %g\n", sinu, sinv, u, v);
      }
    }
  }
}

void calculardxdy (funcL *fL){
  double lobs, mobs, obsra, obsdec, raimage, decimage, ra_offset = 0, dx, dy;
  int i;

  raimage = fL->fg_image->crval[0] * RPDEG;
  decimage = fL->fg_image->crval[1] * RPDEG;

  fL->x0 = (int *) malloc (fL->n_archivos * sizeof (int));
  if (fL->x0 == NULL) {
    fprintf (stderr, "ERROR en calcularSin, fL->x0 = NULL.\n");
    exit (1);
  }
  fL->y0 = (int *) malloc (fL->n_archivos * sizeof (int));
  if (fL->y0 == NULL) {
    fprintf (stderr, "ERROR en calcularSin, fL->y0 = NULL.\n");
    exit (1);
  }
  fL->dx = (double *) malloc (fL->n_archivos * sizeof (double));
  if (fL->dx == NULL) {
    fprintf (stderr, "ERROR en calcularSin, fL->dx = NULL.\n");
    exit (1);
  }
  fL->dy = (double *) malloc (fL->n_archivos * sizeof (double));
  if (fL->dy == NULL) {
    fprintf (stderr, "ERROR en calcularSin, fL->dy = NULL.\n");
    exit (1);
  }

  for (i = 0; i < fL->n_archivos; i++) {
    obsra = ra_offset + fL->header_obs[i]->obsra * RPDEG;
    obsdec = fL->header_obs[i]->obsdec * RPDEG;  
    direccos(obsra, obsdec, raimage, decimage, &lobs, &mobs);
    dx = fL->fg_image->cdelt[0] * RPDEG;   // en radianes
    dy = fL->fg_image->cdelt[1] * RPDEG;   // en radianes
    fL->x0[i] = (fL->fg_image->crpix[0]-1.0) + lobs/dx;
    fL->y0[i] = (fL->fg_image->crpix[1]-1.0) + mobs/dy;
    //L->dx[i] = fL->fg_image->cdelt[0] * 60;   // en arcmin
    //L->dy[i] = fL->fg_image->cdelt[1] * 60;   // en arcmin
    fL->dx[i] = dx;   // en radianes
    fL->dy[i] = dy;   // en radianes
  }
}

void dVdx (funcL *fL, int arch, int iff, int samp, 
	   PoligonoVoronoi *pol, 
	   double *dVdx_R, double *dVdx_I, double *dVdy_R, double *dVdy_I,
	   int aprox) {
  AristaVoronoi *a;
  double integralRx = 0, integralIx = 0, integralRy = 0, integralIy = 0;
  double s0, c1, c2, M, b, u, v, cos_alfa, sen_alfa, x1, x2, y1, y2, xa, ya;
  PoligonoVoronoi *pol_opuesto;
  int nx, ny;
 
  nx = fL->fg_image->size[0];
  ny = fL->fg_image->size[1];
 
  //if (pol->id == 306) {  
  //printf ("dVdx pol %d, arch = %d, iff = %d, samp = %d\n", pol->id, arch, iff, samp);
  //}
  u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
  v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
  (*dVdx_R) = 0;
  (*dVdx_I) = 0;
  (*dVdy_R) = 0;
  (*dVdy_I) = 0;
  //a = (pol->a->poliDer == pol)? pol->a->cwPred : pol->a->cwSucc;
  a = pol->a;

  do {
    if (a->poliDer == pol) {
      cos_alfa = a->cosDer;
      sen_alfa = a->sinDer;
      pol_opuesto = a->poliIzq;
    }
    else if (a->poliIzq == pol) {
      cos_alfa = a->cosIzq;
      sen_alfa = a->sinIzq;
      pol_opuesto = a->poliDer;
    }
    else {
      fprintf (stderr, "ERROR en integraldVdx: arista no coincide con poligono.\n");
      imprimirMallaArchivo (fL->malla, "MallaErrorIntegraldVdx.dat");
      printf ("pol = (%g, %g, %g)\n", pol->x, pol->y, pol->valor);
      imprimirArista (a);
      exit (1);
    }

    xa = (a->ptoIni->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
    ya = (a->ptoIni->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];

    s0 = -xa * cos_alfa  + ya * sen_alfa;
    c1 = -u * cos_alfa + v * sen_alfa;
    c2 = v * cos_alfa + u * sen_alfa;
    if (c2 == 0 || c1 == 0) {
      fprintf (stderr, "u = %g, v = %g, sen_alfa = %g, cos_alfa = %g\n", 
	       u, v, sen_alfa, cos_alfa);
    }
    //s0 = xa * cos_alfa  + ya * sen_alfa;
    //c1 = u * cos_alfa + v * sen_alfa;
    //c2 = v * cos_alfa - u * sen_alfa;
    x1 = (pol->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
    y1 = (pol->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
    x2 = (pol_opuesto->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
    y2 = (pol_opuesto->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];

    // Para dVdx
    if (x2 != x1) {
      M = - cos_alfa * sen_alfa / (x2 - x1);
      b = cos_alfa * (s0 * cos_alfa + x1) / (x2 - x1);   
    }
    else {
      M = sen_alfa * sen_alfa / (y2 - y1);
      b = - sen_alfa * (s0 * cos_alfa + x1) / (y2 - y1);   
    }
    /*
    if (x2 != x1) {
      M = - cos_alfa * sen_alfa / (x2 - x1);
      b = cos_alfa * (s0 * cos_alfa - x1) / (x2 - x1);   
    }
    else {
      M = - sen_alfa * sen_alfa / (y2 - y1);
      b = sen_alfa * (s0 * cos_alfa - x1) / (y2 - y1);   
    }
    */
    if (aprox == APROX_DX) {
      integraldVdxAprox (fL, arch, iff, samp, pol, a, 
			 s0, c1, c2, M, b, &integralRx, &integralIx);
    }
    else {
      integraldVdx (fL, arch, iff, samp, pol, a, 
		    s0, c1, c2, M, b, &integralRx, &integralIx);
    }
    // Para dVdy
    if (y2 != y1) {
      M = cos_alfa * sen_alfa / (y2 - y1);
      b = sen_alfa * (s0 * sen_alfa - y1) / (y2 - y1);   
    }
    else {
      M = cos_alfa * cos_alfa / (x1 - x2);
      b = cos_alfa * (s0 * sen_alfa - y1) / (x1 - x2);   
    }
    /*
    if (y2 != y1) {
      M = cos_alfa * sen_alfa / (y2 - y1);
      b = sen_alfa * (s0 * sen_alfa - y1) / (y2 - y1);   
    }
    else {
      M = cos_alfa * cos_alfa / (x2 - x1);
      b = cos_alfa * (s0 * sen_alfa - y1) / (x2 - x1);   
    }
    */
    if (aprox == APROX_DX) {
      integraldVdxAprox (fL, arch, iff, samp, pol, a, 
			 s0, c1, c2, M, b, &integralRy, &integralIy);
    }
    else {
      integraldVdx (fL, arch, iff, samp, pol, a, 
		    s0, c1, c2, M, b, &integralRy, &integralIy);
    }
    //if (pol->id == 306) {
    //printf ("ok\n");
    //}
    //integralR *= fL->fg_scale;
    //integralI *= fL->fg_scale;
    //(*dVdx_R) -= (pol->valor - pol_opuesto->valor) * integralRx;
    (*dVdx_R) += (pol->valor - pol_opuesto->valor) * integralRx;
    (*dVdy_R) += (pol->valor - pol_opuesto->valor) * integralRy;
    //(*dVdx_I) -= (pol->valor - pol_opuesto->valor) * integralIx;
    (*dVdx_I) += (pol->valor - pol_opuesto->valor) * integralIx;
    (*dVdy_I) += (pol->valor - pol_opuesto->valor) * integralIy;

    if (a->poliDer == pol) {
      a = a->cwPred;
    }
    else if (a->poliIzq == pol) {
      a = a->cwSucc;
    }
    else {
      fprintf (stderr, "ERROR en dVdx: Mal recorrido de aristas.\n");
      fprintf (stderr, "Arista:\n");
      imprimirArista (a);
      imprimirMallaArchivo(fL->malla, "MallaErrordVdx.dat");
      if (a->poliDer == NULL) {
	fprintf (stderr, "a->poliDer = NULL\n");
      }
      else {
	fprintf (stderr, "a->poliDer = (%g, %g)\n", 
		 a->poliDer->x, a->poliDer->y);
      }      
      if (a->poliIzq == NULL) {
	fprintf (stderr, "a->poliIzq = NULL\n");
      }
      else {
	fprintf (stderr, "a->poliIzq = (%g, %g)\n", 
		 a->poliIzq->x, a->poliIzq->y);
      }
      fprintf (stderr, "pol = (%g, %g)\n", pol->x, pol->y);
      exit (1);
    }
  }
  while (a != pol->a);
  /*
  if ((*dVdx_R) == 0 || (*dVdx_I) == 0 || (*dVdy_R) == 0 || (*dVdy_I) == 0) {
    fprintf (stderr, "ERROR: dVdx = (%g, %g), dVdy = (%g, %g)\n", 
	     (*dVdx_R), (*dVdx_I), (*dVdy_R), (*dVdy_I));
    a = pol->a;
    do {
      integraldVdx (fL, arch, iff, samp, pol, a, &integralR, &integralI);
      fprintf (stderr, "integralR = %g, integralI = %g\n", 
	       integralR, integralI);
      if (a->poliDer == pol) {
	a = a->cwPred;
      }
      else if (a->poliIzq == pol) {
	a = a->cwSucc;
      }
    } while (a != pol->a);
    exit (0);
  }
  */
  (*dVdx_R) *= fL->fg_scale * fabs (nx * fL->dx[arch]);
  (*dVdx_I) *= fL->fg_scale * fabs (nx * fL->dx[arch]);
  (*dVdy_R) *= fL->fg_scale * fabs (ny * fL->dy[arch]);
  (*dVdy_I) *= fL->fg_scale * fabs (ny * fL->dy[arch]);
  //printf ("dVdx_I = %g\n", *dVdx_I);
}

void integraldVdx (funcL *fL, int arch, int iff, int samp, 
		   PoligonoVoronoi *pol, AristaVoronoi *a,
		   double s0, double c1, double c2, double M, double b,
		   double *integralR, double *integralI) {
  int i, j, nx, ny, tipoLinea, lineaAnterior = 0;
  PuntoVoronoi *p1, *p2, *pIni, *pFin, *paux;
  double m, u, v, I1R, I1I;
  double x1, y1, x2, y2, sin_alfa, cos_alfa, t, dt, sen1, cos1, sen2, cos2;
  
  p1 = newPuntoVoronoi (0, 0);
  p2 = newPuntoVoronoi (0, 0);
  nx = fL->fg_image->size[0];
  ny = fL->fg_image->size[1];
  (*integralR) = 0;
  (*integralI) = 0;
  u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
  v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
  //printf ("1) u, v = (%g, %g)\n", 
  //fL->samples_mod[arch][samp].u, fL->samples_mod[arch][samp].v);
  
  m = (a->ptoIni->y - a->ptoFin->y) * ny / 
    ((a->ptoIni->x - a->ptoFin->x) * nx);
  if (m > 1) {
    tipoLinea = 3;
  }
  else if (m > 0) {
    tipoLinea = 1;
  }
  else if (m > -1) {
    tipoLinea = 2;
  }
  else {
    tipoLinea = 4;
  }

  if (a->ptoIni->x < a->ptoFin->x) {
    if (!interseccionCuadrado (a->ptoIni, a->ptoFin, p1, p2, tipoLinea)) {
      free (p1);
      free (p2);
      //printf ("no topa\n");
      return;
    }
  }
  else {
    if (!interseccionCuadrado (a->ptoFin, a->ptoIni, p2, p1, tipoLinea)) {
      free (p1);
      free (p2);
      //printf ("no topa\n");
      return;
    }
  }
  
  //printf ("p1 = (%g, %g), p2 = (%g, %g)\n", p1->x, p1->y, p2->x, p2->y);

  if (a->poliDer == pol) {
    paux = p1;
    p1 = p2;
    p2 = paux;
    cos_alfa = a->cosDer;
    sin_alfa = a->sinDer;
  }
  else if (a->poliIzq == pol) {
    cos_alfa = a->cosIzq;
    sin_alfa = a->sinIzq;
  }
  else {
    fprintf (stderr, "ERROR en integraldVdx: arista no coincide con poligono.\n");
    imprimirMallaArchivo (fL->malla, "MallaErrorIntegraldVdx.dat");
    printf ("pol = (%g, %g, %g)\n", pol->x, pol->y, pol->valor);
    imprimirArista (a);
    if (a->poliDer == NULL) {
      printf ("poliDer = NULL\n");
    }
    else {
      printf ("poliDer = (%g, %g, %g)\n", 
	      a->poliDer->x, a->poliDer->y, a->poliDer->valor);
    }
    if (a->poliIzq == NULL) {
      printf ("poliIzq = NULL\n");
    }
    else {
      printf ("poliIzq = (%g, %g, %g)\n", 
	      a->poliIzq->x, a->poliIzq->y, a->poliIzq->valor);
    }
    exit (1);
  }

  //printf ("p1 = (%g, %g), p2 = (%g, %g)\n", p1->x, p1->y, p2->x, p2->y);
  //printf ("cosalfa = %g, sinalfa = %g\n", cos_alfa, sin_alfa);

  // Damos vuelta la imagen
  /* 
  paux = p1;
  p1 = p2;
  p2 = paux;
  cos_alfa = -cos_alfa;
  */

  // Integramos de abajo hacia arriba
  /*
  if (p1->y > p2->y) {
    paux = p1;
    p1 = p2;
    p2 = paux;
    //printf ("abajo hacia arriba\n");
  }
  */
  //p1->x = 1 - p1->x;
  //p2->x = 1 - p2->x;

  //x1 = (pIni->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
  //y1 = (pIni->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
  x1 = (p1->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
  y1 = (p1->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
  //s0 = y1 * sin_alfa + x1 * cos_alfa;
  //c1 = PI * (v * cos_alfa - u * sin_alfa);
  //c2 = PI * s0 * (u * cos_alfa + v * sin_alfa);

  pFin = newPuntoVoronoi (p1->x, p1->y);
  
  for (pIni = newPuntoVoronoi (p1->x, p1->y);
       (pFin->x != p2->x) || (pFin->y != p2->y);) {
    //int linant = lineaAnterior;
    lineaAnterior = encontrarSgtePtoPixel (fL, pIni, p2, pFin, lineaAnterior);
    /*
    if (fabs (pIni->x - pFin->x) < 1e-6 && 
	fabs (pIni->y - pFin->y) < 1e-6) {
      fprintf (stderr, "ERROR integraldVdx: pIni = Pfin\n");
      fprintf (stderr, "      pIni = (%g, %g), pFin = (%g, %g)\n", 
	       pIni->x, pIni->y, pFin->x, pFin->y);
      fprintf (stderr, "      errorx = %.10g, errory = %.10g\n", 
	       pIni->x - pFin->x, pIni->y - pFin->y);
      exit (1);
    }
    */
    if (SQR(pIni->x - pFin->x) + SQR(pIni->y - pFin->y) > 
	SQR(1.0/(nx - 1.0)) + SQR(1.0/(ny - 1.0)) + 1e-6) {
      fprintf (stderr, "ERROR integraldVdx: d(pIni, pFin) muy grande.\n");
      fprintf (stderr, "      pIni = (%g, %g), pFin = (%g, %g)\n", 
	       pIni->x, pIni->y, pFin->x, pFin->y);
      fprintf (stderr, "      d = %.10g vs %.10g\n", 
	       SQR(pIni->x - pFin->x) + SQR(pIni->y - pFin->y), 
	       SQR(1.0/(nx - 1.0)) + SQR(1.0/(ny - 1.0)));
      fprintf (stderr, "      dx = %.10g, dy =  %.10g (en pixeles)\n", 
	       (pFin->x - pIni->x) * (nx - 1), (pFin->y - pIni->y) * (ny - 1));
      fprintf (stderr, "      pIni = (%g, %g) pFin =  (%g, %g) (en pixeles)\n",
	       pIni->x * (nx - 1), pIni->y * (nx - 1), 
	       pFin->x * (ny - 1), pFin->y * (ny - 1));
      fprintf (stderr, "      lineaAnterior = %d\n", lineaAnterior);
      
      exit (1);
      }
    
    x1 = (pIni->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
    y1 = (pIni->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
    x2 = (pFin->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
    y2 = (pFin->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
    
    t = 0.5 * (cos_alfa * (y1 + y2) + sin_alfa * (x1 + x2));
    dt = sin_alfa * (x2 - x1) + cos_alfa * (y2 - y1);
    //t = 0.5 * (cos_alfa * (y1 + y2) - sin_alfa * (x1 + x2));
    //dt = sin_alfa * (x1 - x2) + cos_alfa * (y2 - y1);
    cos1 = cos (2 * PI * (c1*s0 + c2*t));
    sen1 = sin (2 * PI * (c1*s0 + c2*t));
    cos2 = cos (PI * c2 * dt);
    sen2 = sin (PI * c2 * dt);
    //printf (" sin_c1dt = %g, c1 = %g\n", sin_c1dt, c1);
    // Trans Fourier

    if (c2 != 0) {
      I1R = (cos1 * (M * t + b) * sen2 - 
	     sen1 * M * 0.5 * (sen2 /(PI * c2) - dt * cos2)) / (PI * c2);
      I1I = (sen1 * (M * t + b) * sen2 + 
	     cos1 * M * 0.5 * (sen2 /(PI * c2) - dt * cos2)) / (PI * c2);
    }
    else {
      I1R = (cos1 * (M * t + b) * dt);
      I1I = (sen1 * (M * t + b) * dt);
    }
    // Trans Inversa
    /*
    if (c1 != 0) {
      I1R = cos (-2 * (c2 + c1 * t)) * sin_c1dt / c1;
      I1I = sin (-2 * (c2 + c1 * t)) * sin_c1dt / c1;
    } else {
      I1R = cos (-2 * (c2 + c1 * t));
      I1I = sin (-2 * (c2 + c1 * t));
    }
    */
    i = round (0.5 * (pIni->x + pFin->x) * (nx - 1.0));
    j = round (0.5 * (pIni->y + pFin->y) * (ny - 1.0));
    (*integralR) += fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff] * I1R;
    (*integralI) += fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff] * I1I;

    if (isnan((*integralR)) || isnan((*integralI))) {
      fprintf (stderr, "ERROR en integraldVdx\n");
      fprintf (stderr, "      cos1 = %g, sen1 = %g, cos2 = %g, sen2 = %g\n", 
	       cos1, sen1, cos2, sen2);
      fprintf (stderr, "M = %g, t = %g, b = %g, dt = %g, c2 = %g\n", M, t, b, dt, c2);
      fprintf (stderr, "      I1R = %g, I1I = %g, atten = %g\n", 
	       I1R, I1I, fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff]);
      fprintf (stderr, "      arch = %d, i = %d, j = %d size = %d, iff = %d\n", 
	       arch, i, j, fL->fg_image->size[0], iff);
      return;
    }
    fL->PixIntegrados_x[i][j]++;
    
    pIni->x = pFin->x;
    pIni->y = pFin->y;
  }

  //(*integralR) /= (PI * c2);
  //(*integralI) /= (PI * c2);
  /*
  if ((*integralR) != 0 && (*integralI) != 0) {
    t = 0.5 * (cos_alfa * (p1->y + p2->y) - sin_alfa * (p1->x + p2->x));
    dt = sin_alfa * (p1->x - p2->x) + cos_alfa * (p2->y - p1->y);
    printf ("s0 = %g, c1 = %g, c2 = %g, M = %g, b= %g, t = %g, dt = %g\n", 
	    s0, c1, c2, M, b, t, dt);
    printf ("IR = %g, II = %g\n------------\n", (*integralR), (*integralI));
  }
  */
  //(*integralR) *= -sin_alfa;
  //(*integralI) *= -sin_alfa;
  free (pIni);
  free (pFin);
  
  free (p1);
  free (p2);
}

void integraldVdxAprox (funcL *fL, int arch, int iff, int samp, 
			PoligonoVoronoi *pol, AristaVoronoi *a,
			double s0, double c1, double c2, double M, double b,
			double *integralR, double *integralI) {
  int i, j, nx, ny, tipoLinea, lineaAnterior = 0;
  PuntoVoronoi *p1, *p2, *pIni, *pFin, *paux;
  double m, I1R, I1I;
  double x1, y1, x2, y2, t, dt, sin_alfa, cos_alfa;
  
  p1 = newPuntoVoronoi (0, 0);
  p2 = newPuntoVoronoi (0, 0);
  nx = fL->fg_image->size[0];
  ny = fL->fg_image->size[1];
  (*integralR) = 0;
  (*integralI) = 0;
  //u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
  //v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
  //printf ("1) u, v = (%g, %g)\n", 
  //fL->samples_mod[arch][samp].u, fL->samples_mod[arch][samp].v);
  //printf ("2) u, v = (%g, %g)\n", 
	  //u, v);
  
  m = (a->ptoIni->y - a->ptoFin->y) * ny / 
    ((a->ptoIni->x - a->ptoFin->x) * nx);
  if (m > 1) {
    tipoLinea = 3;
  }
  else if (m > 0) {
    tipoLinea = 1;
  }
  else if (m > -1) {
    tipoLinea = 2;
  }
  else {
    tipoLinea = 4;
  }

  if (a->ptoIni->x < a->ptoFin->x) {
    if (!interseccionCuadrado (a->ptoIni, a->ptoFin, p1, p2, tipoLinea)) {
      free (p1);
      free (p2);
      //printf ("no topa\n");
      return;
    }
  }
  else {
    if (!interseccionCuadrado (a->ptoFin, a->ptoIni, p2, p1, tipoLinea)) {
      free (p1);
      free (p2);
      //printf ("no topa\n");
      return;
    }
  }

  if (a->poliDer == pol) {
    paux = p1;
    p1 = p2;
    p2 = paux;
    cos_alfa = a->cosDer;
    sin_alfa = a->sinDer;
  }
  else if (a->poliIzq == pol) {
    cos_alfa = a->cosIzq;
    sin_alfa = a->sinIzq;
  }
  else {
    fprintf (stderr, "ERROR en integraldVdx: arista no coincide con poligono.\n");
    imprimirMallaArchivo (fL->malla, "MallaErrorIntegraldVdx.dat");
    printf ("pol = (%g, %g, %g)\n", pol->x, pol->y, pol->valor);
    imprimirArista (a);
    if (a->poliDer == NULL) {
      printf ("poliDer = NULL\n");
    }
    else {
      printf ("poliDer = (%g, %g, %g)\n", 
	      a->poliDer->x, a->poliDer->y, a->poliDer->valor);
    }
    if (a->poliIzq == NULL) {
      printf ("poliIzq = NULL\n");
    }
    else {
      printf ("poliIzq = (%g, %g, %g)\n", 
	      a->poliIzq->x, a->poliIzq->y, a->poliIzq->valor);
    }
    exit (1);
  }

  x1 = (p1->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
  y1 = (p1->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];

  pFin = newPuntoVoronoi (p1->x, p1->y);

  /*
  if (p1->x < p2->x) {
    signo = 1;
  }
  else {
    signo = -1;
  }
  */

  for (pIni = newPuntoVoronoi (p1->x, p1->y);
       (pFin->x != p2->x) || (pFin->y != p2->y);) {
    //int linant = lineaAnterior;

    lineaAnterior = encontrarSgtePtoPixel (fL, pIni, p2, pFin, lineaAnterior);
    //printf ("  sumando {(%g, %g), (%g, %g)}\n", 
    //pIni->x * (nx - 1), pIni->y * (ny - 1), 
    //pFin->x * (nx - 1), pFin->y * (ny - 1));

    if (SQR(pIni->x - pFin->x) + SQR(pIni->y - pFin->y) > 
    SQR(1.0/(nx - 1.0)) + SQR(1.0/(ny - 1.0)) + 1e-6) {
    //if (pol->id == 306) {
      fprintf (stderr, "ERROR integraldVdx: d(pIni, pFin) muy grande.\n");
      fprintf (stderr, "      pIni = (%g, %g), pFin = (%g, %g)\n", 
	       pIni->x, pIni->y, pFin->x, pFin->y);
      fprintf (stderr, "      d = %.10g vs %.10g\n", 
	       SQR(pIni->x - pFin->x) + SQR(pIni->y - pFin->y), 
	       SQR(1.0/(nx - 1.0)) + SQR(1.0/(ny - 1.0)));
      fprintf (stderr, "      dx = %.10g, dy =  %.10g (en pixeles)\n", 
	       (pFin->x - pIni->x) * (nx - 1), (pFin->y - pIni->y) * (ny - 1));
      fprintf (stderr, "      pIni = (%g, %g) pFin =  (%g, %g) (en pixeles)\n",
	       pIni->x * (nx - 1), pIni->y * (nx - 1), 
	       pFin->x * (ny - 1), pFin->y * (ny - 1));
      fprintf (stderr, "      lineaAnterior = %d\n", lineaAnterior);
      
      exit (1);
      }
    
    x1 = (pIni->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
    y1 = (pIni->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
    x2 = (pFin->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
    y2 = (pFin->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];

    t = 0.5 * (cos_alfa * (y1 + y2) + sin_alfa * (x1 + x2));
    dt = sin_alfa * (x2 - x1) + cos_alfa * (y2 - y1);
    //t = 0.5 * (cos_alfa * (y1 + y2) - sin_alfa * (x1 + x2));
    //dt = sin_alfa * (x1 - x2) + cos_alfa * (y2 - y1);
    /*
    if (dt < -0) {
      printf ("ERROR integralAprox: dt = %g\n", dt);
      printf ("                     x1 = %g, x2= %g\n", x1, x2);
      printf ("                     y1 = %g, y2= %g\n", y1, y2);
      printf ("                     sin = %g, cos = %g\n\n", sin_alfa, cos_alfa);
      //exit (1);
    }
    */

    I1R = dt * (M*t + b) * cos (2 * PI * (t * c2 + s0 * c1));
    I1I = dt * (M*t + b) * sin (2 * PI * (t * c2 + s0 * c1));
  
    // Trans Inversa
    /*
      I1R = cos (-2 * PI *(u * x + v * y)) * d;
      I1I = sin (-2 * PI *(u * x + v * y)) * d;
    */

    i = round (0.5 * (pIni->x + pFin->x) * (nx - 1.0));
    j = round (0.5 * (pIni->y + pFin->y) * (ny - 1.0));
    (*integralR) += fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff] * I1R;
    (*integralI) += fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff] * I1I;
    fL->PixIntegrados_x[i][j]++;
    
    pIni->x = pFin->x;
    pIni->y = pFin->y;
  }

  //(*integralR) *= -sin_alfa;
  //(*integralI) *= -sin_alfa;
  free (pIni);
  free (pFin);
  
  free (p1);
  free (p2);
}


/* Encuentra la coordenada de la siguiente interseccion de la arista
 * con el pixel, donde p1 es el punto que intersecta con el borde
 * anterior, pFin el ultimo punto de la arista, pSgte el punto a
 * encontrar y linea anterior es el borde donde topo el punto
 * anterior.
 */
int encontrarSgtePtoPixel (funcL *fL, PuntoVoronoi *p1, PuntoVoronoi *pFin,
			   PuntoVoronoi *pSgte, int lineaAnterior) {
  double x, y, dx1, dy1, d, d1, d2, error2, delta = 1e-6;
  int i, j, nx, ny, ret = 0;
  
  nx = fL->fg_image->size[0];
  ny = fL->fg_image->size[1];
  //error2 = pow (1e-3 / (nx - 1.0), 2) + pow (1e-3 / (ny - 1.0), 2); 
  error2 = SQR(1e-3 / (nx - 1.0)) + SQR (1e-3 / (ny - 1.0)); 
  dx1 = pFin->x - p1->x;
  dy1 = pFin->y - p1->y;
  d = 1e300;

  // para evitar problemas numericos...
  switch (lineaAnterior) {
  case BORDE_IZQ:
    i = round (p1->x * (nx - 1) - delta);
    j = round (p1->y * (ny - 1));
    //printf ("caso BORDE_IZQ\n");
    break;
  case BORDE_DER:
    i = round (p1->x * (nx - 1) + delta);
    j = round (p1->y * (ny - 1));
    //printf ("caso BORDE_DER\n");
    break;
  case BORDE_INF:
    i = round (p1->x * (nx - 1));
    j = round (p1->y * (ny - 1) - delta);
    //printf ("caso BORDE_INF\n");
    break;
  case BORDE_SUP:
    i = round (p1->x * (nx - 1));
    j = round (p1->y * (ny - 1) + delta);
    //printf ("caso BORDE_SUP\n");
    break;
  default:
    i = round (p1->x * (nx - 1));
    j = round (p1->y * (ny - 1));
    //printf ("caso 0\n");
    break;
  }
  /*
  if (i >= 76 && i <= 77 && j >= 25 && j <= 26) {
  printf ("lineaAnterior = %d, (i, j) = (%d, %d)\n", lineaAnterior, i, j);
  printf ("p1 = (%g, %g), pFin = (%g, %g)\n", p1->x, p1->y, pFin->x, pFin->y);
  printf ("p1 = (%g, %g), pFin = (%g, %g) (en pixeles)\n", 
	  p1->x * (nx - 1), p1->y * (ny - 1), 
	  pFin->x * (nx - 1), pFin->y * (ny - 1));
  printf ("p1 - (i, j) = (%g, %g)\n", 
	  p1->x * (nx - 1), p1->y * (ny - 1) - j - 0.5);
  }
  */

  // linea inferior.
  if (lineaAnterior != BORDE_SUP) {
    y = (j - 0.5) / (ny - 1.0);
    x = (y - p1->y) * dx1 / dy1 + p1->x;
    //printf ("inferior y = %g, %g en pixeles, ppunto = %g\n", y, y * (ny - 1),
    //(x - p1->x) * (pFin->x - p1->x) + 
    //(y - p1->y) * (pFin->y - p1->y));
    if ((x - p1->x) * (pFin->x - p1->x) + 
	(y - p1->y) * (pFin->y - p1->y) >= 0) {
      //(y - p1->y) * (pFin->y - p1->y) > 0 - delta / (nx - 1.0)) {
      //d1 = pow (p1->x - x, 2) + pow (p1->y - y, 2);
      d1 = SQR (p1->x - x) + SQR (p1->y - y);
      //printf ("OK1 d1 = %g, d = %g\n", d1, d);
      //printf ("inf d1 = %g, d = %g\n", d1, d);
      if (d1 * nx * ny < SQR (delta)) {
	y = (j - 1.5) / (ny - 1.0);
	x = (y - p1->y) * dx1 / dy1 + p1->x;
	if ((x - p1->x) * (pFin->x - p1->x) + 
	    (y - p1->y) * (pFin->y - p1->y) >= 0) {
	  //(y - p1->y) * (pFin->y - p1->y) > 0 - delta / (nx - 1.0)) {
	  d1 = SQR (p1->x - x) + SQR (p1->y - y);
	  //printf ("inf d1 = %g, d = %g\n", d1, d);
	}
      }
      if (d1 < d && d1 * nx * ny > SQR (delta)) {
	d = d1;
	pSgte->x = x;
	pSgte->y = y;
	ret = BORDE_INF;
	//printf ("OK, (%g, %g) (en pixeles)\n", x * (nx - 1), y * (ny - 1));
      }
    }
  }
  
  // linea derecha.
  if (lineaAnterior != BORDE_IZQ) {
    x = (i + 0.5) / (nx - 1.0);
    y = (x - p1->x) * dy1 / dx1 + p1->y;
    //printf ("derecha x = %g, ppunto = %g\n", x, 
    //    (x - p1->x) * (pFin->x - p1->x) + 
    //    (y - p1->y) * (pFin->y - p1->y));
    if ((x - p1->x) * (pFin->x - p1->x) + 
	(y - p1->y) * (pFin->y - p1->y) >= 0) {
      //(y - p1->y) * (pFin->y - p1->y) > 0 - delta / (nx - 1.0)) {
      //d1 = pow (p1->x - x, 2) + pow (p1->y - y, 2); 
      d1 = SQR (p1->x - x) + SQR (p1->y - y); 
      //printf ("OK1 d1 = %g, d = %g\n", d1, d);
      //printf ("der d1 = %g, d = %g\n", d1, d);
      if (d1 * nx * ny < SQR (delta)) {
	x = (i + 1.5) / (nx - 1.0);
	y = (x - p1->x) * dy1 / dx1 + p1->y;
	if ((x - p1->x) * (pFin->x - p1->x) + 
	    (y - p1->y) * (pFin->y - p1->y) >= 0) {
	  //(y - p1->y) * (pFin->y - p1->y) > 0 - delta / (nx - 1.0)) {
	  d1 = SQR (p1->x - x) + SQR (p1->y - y);
	  //printf ("der d1 = %g, d = %g\n", d1, d);
	}
      }
      if (d1 < d && d1 * nx * ny > SQR (delta)) {
	d = d1;
	pSgte->x = x;
	pSgte->y = y;
	ret = BORDE_DER;
	//printf ("OK\n");
      }
    }
  }

  // linea superior.
  if (lineaAnterior != BORDE_INF) {
    y = (j + 0.5) / (ny - 1.0);
    x = (y - p1->y) * dx1 / dy1 + p1->x;
    //printf ("superior y = %g, ppunto = %g\n", y, 
    //    (x - p1->x) * (pFin->x - p1->x) + 
    //    (y - p1->y) * (pFin->y - p1->y));
    if ((x - p1->x) * (pFin->x - p1->x) + 
	(y - p1->y) * (pFin->y - p1->y) >= 0) {
      //(y - p1->y) * (pFin->y - p1->y) > 0 - delta / (nx - 1.0)) {
      //d1 = pow (p1->x - x, 2) + pow (p1->y - y, 2); 
      d1 = SQR (p1->x - x) + SQR (p1->y - y); 
      //printf ("OK1 d1 = %g, d = %g\n", d1, d);
      //printf ("sup d1 = %g, d = %g\n", d1, d);
      if (d1 * nx * ny < SQR (delta)) {
	y = (j + 1.5) / (ny - 1.0);
	x = (y - p1->y) * dx1 / dy1 + p1->x;
	if ((x - p1->x) * (pFin->x - p1->x) + 
	    (y - p1->y) * (pFin->y - p1->y) >= 0) {
	  //(y - p1->y) * (pFin->y - p1->y) > 0 - delta / (nx - 1.0)) {
	  d1 = SQR (p1->x - x) + SQR (p1->y - y);
	  //printf ("sup d1 = %g, d = %g\n", d1, d);
	}
      }
      if (d1 < d && d1 * nx * ny > SQR (delta)) {
	d = d1;
	pSgte->x = x;
	pSgte->y = y;
	ret = BORDE_SUP;
	//printf ("OK\n");
      }
    }
  }

  // linea izquierda.
  if (lineaAnterior != BORDE_DER) {
    x = (i - 0.5) / (nx - 1.0);
    y = (x - p1->x) * dy1 / dx1 + p1->y;
    //printf ("izquierda x = %g, ppunto = %.10g\n", x, 
    //(x - p1->x) * (pFin->x - p1->x) + 
    //(y - p1->y) * (pFin->y - p1->y));
    //printf ("(x, y) = (%g, %g)\n", x * (nx - 1), y * (ny - 1));
    if ((x - p1->x) * (pFin->x - p1->x) + 
	(y - p1->y) * (pFin->y - p1->y) >= 0) {
      //(y - p1->y) * (pFin->y - p1->y) > 0 - delta / (nx - 1.0)) {
      //d1 = pow (p1->x - x, 2) + pow (p1->y - y, 2); 
      d1 = SQR (p1->x - x) + SQR (p1->y - y); 
      //printf ("OK1 d1 = %g, d = %g\n", d1, d);
      //printf ("izq d1 = %g, d = %g\n", d1, d);
      if (d1 * nx * ny < SQR (delta)) {
	x = (i - 1.5) / (nx - 1.0);
	y = (x - p1->x) * dy1 / dx1 + p1->y;
	if ((x - p1->x) * (pFin->x - p1->x) + 
	    (y - p1->y) * (pFin->y - p1->y) >= 0) {
	  //(y - p1->y) * (pFin->y - p1->y) > 0 - delta / (nx - 1.0)) {
	  d1 = SQR (p1->x - x) + SQR (p1->y - y);
	  //printf ("izq d1 = %g, d = %g\n", d1, d);
	}
      }
      if (d1 < d && d1 * nx * ny > SQR (delta)) {
	d = d1;
	pSgte->x = x;
	pSgte->y = y;
	ret = BORDE_IZQ;
	//printf ("OK\n");
      }
    }
  }

  // comparacion con pFin.
  //d2 = pow(p1->x - pFin->x, 2) + pow(p1->y - pFin->y, 2);
  d2 = SQR(p1->x - pFin->x) + SQR(p1->y - pFin->y);
  if (d >= d2 && d1 * nx * ny > SQR (delta)) {
    pSgte->x = pFin->x;
    pSgte->y = pFin->y;
    //printf ("OK Fin, d = %g, d2 = %g\n", d, d2);
  }
  //printf ("p1 = (%g, %g), pSgte = (%g, %g)\n", 
  //p1->x * (nx - 1), p1->y * (ny - 1), 
  //pSgte->x * (nx - 1), pSgte->y * (ny - 1));
  //printf ("ultimo OK, d = %g\n\n", d);
  return ret;
}

void calcularFourierI (funcL *fL) {
  int i, j, arch, iff, samp;
  double Re, Im, x, y, kx, u, v;
  
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      fL->fourierI[i][j] = 0;
      fL->nPixIntegrados_I++; 
      fL->PixIntegrados_I[i][j]++;

      for (arch = 0; arch < fL->n_archivos; arch++) {
	y = (j - fL->y0[arch]) * fL->dy[arch];
	x = (i - fL->x0[arch]) * fL->dx[arch];
	//printf ("nstokes = %d\n", fL->header_obs[arch]->nstokes);
	for (iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
	  double temp = 0;
	  for(samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	    u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
	    v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
	    kx = (u * x + v * y);
	    /*
	    printf ("VISobs %g + i %g\n", 
		    fL->samples_obs[arch][samp].rdata[3 * iff],
		    fL->samples_obs[arch][samp].rdata[3 * iff + 1]);
	    printf ("     w %g\n",
		    fL->samples_obs[arch][samp].rdata[3 * iff + 2]);
	    printf ("   mod %g + i %g\n",
		    fL->samples_mod[arch][samp].rdata[3 * iff],
		    fL->samples_mod[arch][samp].rdata[3 * iff + 1]);
	    printf ("     w %g\n",
		    fL->samples_mod[arch][samp].rdata[3 * iff + 2]);
	    */
	    // Transf Fourier Directa

	    Re = ((fL->samples_mod[arch][samp].rdata[3 * iff] 
		   - fL->samples_obs[arch][samp].rdata[3 * iff]) 
		  * cos (2 * PI * kx));
	    Im = ((fL->samples_mod[arch][samp].rdata[3 * iff + 1] 
		   - fL->samples_obs[arch][samp].rdata[3 * iff + 1])
		  * sin (2 * PI * kx));
	    /*
	    // Transf Fourier Inversa
	    /*
	    Re = ((fL->samples_mod[arch][samp].rdata[3 * iff] 
		   - fL->samples_obs[arch][samp].rdata[3 * iff]) 
		  * cos (-2 * PI * kx));
	    Im = ((fL->samples_mod[arch][samp].rdata[3 * iff + 1] 
		   - fL->samples_obs[arch][samp].rdata[3 * iff + 1])
		  * sin (-2 * PI * kx));

	    */
	    //temp += fL->sin[arch][samp][iff] // 10 Mayo 2007
	    //* fL->samples_mod[arch][samp].rdata[3 * iff + 2] * (Re + Im);
	    //fL->fourierI[i][j] += fL->sin[arch][samp][iff] // 29 Oct 2005
	    //* fL->samples_mod[arch][samp].rdata[3 * iff + 2] * (Re + Im);
	    fL->fourierI[i][j] += fL->sin[arch][samp][iff] 
	      * fL->samples_mod[arch][samp].rdata[3 * iff + 2] * (Re + Im)
	      * fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff];
	  }
	  //temp *= fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff]; // 10 Mayo 2007
	  //fL->fourierI[i][j] += temp; // 10 Mayo 2007
	  //fL->fourierI[i][j] *= fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff]; // 29 Oct 2005
	}
      }
    }
  }
}

void calcularFourierIAprox (funcL *fL) {
  int i, j, arch, iff, samp;
  double Re, Im, x, y, kx, u, v, dx, dy;
  
  dx = fL->fg_image->cdelt[0] * RPDEG;   // en radianes
  dy = fL->fg_image->cdelt[1] * RPDEG;   // en radianes
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      fL->fourierI[i][j] = 0;
      fL->nPixIntegrados_I++; 
      fL->PixIntegrados_I[i][j]++;

      for (arch = 0; arch < fL->n_archivos; arch++) {
	y = (j - fL->y0[arch]) * fL->dy[arch];
	x = (i - fL->x0[arch]) * fL->dx[arch];
	for (iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
	  for(samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	    u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
	    v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
	    kx = (u * x + v * y);
	    // Transf Fourier Directa

	    Re = ((fL->samples_mod[arch][samp].rdata[3 * iff] 
		   - fL->samples_obs[arch][samp].rdata[3 * iff]) 
		  * cos (2 * PI * kx));
	    Im = ((fL->samples_mod[arch][samp].rdata[3 * iff + 1] 
		   - fL->samples_obs[arch][samp].rdata[3 * iff + 1])
		  * sin (2 * PI * kx));
	    
	    // Transf Fourier Inversa
	    /*
	    Re = ((fL->samples_mod[arch][samp].rdata[3 * iff] 
		   - fL->samples_obs[arch][samp].rdata[3 * iff]) 
		  * cos (-2 * PI * kx));
	    Im = ((fL->samples_mod[arch][samp].rdata[3 * iff + 1] 
		   - fL->samples_obs[arch][samp].rdata[3 * iff + 1])
		  * sin (-2 * PI * kx));

	    */
	    fL->fourierI[i][j] += fL->samples_mod[arch][samp].rdata[3 * iff + 2] 
	    * (Re + Im);
	    //fL->fourierI[i][j] += fL->samples_mod[arch][samp].rdata[3 * iff + 2] 
	    //* (Re + Im) * fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff];
	  }
	  fL->fourierI[i][j] *= fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff];
	}
      }
      fL->fourierI [i][j] *= fabs(dx * dy);
    }
  }
}

void eliminarFuncL (funcL *fL){
  int i, j, arch, samp;

  for (arch = 0; arch < fL->n_archivos; arch++) {
    for(samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
      free (fL->sin[arch][samp]);
    }
    free (fL->sin[arch]);
  }
  free (fL->sin);
  for (i = 0; i < fL->n_archivos; i++) {
    free(fL->samples_obs[i]);
    free(fL->header_obs[i]);
    free(fL->samples_mod[i]);
    free (fL->infile[i]);
    fL->infile[i] = NULL;
  }
  free (fL->header_obs);
  free (fL->samples_obs);
  free (fL->samples_mod);
  free (fL->infile);
  free (fL->nombreFits);
  free (fL->x0);
  free (fL->y0);
  free (fL->dx);
  free (fL->dy);

  //fL->atten= atenuacion(fL->fg_image, fL->header_obs, fL->n_archivos, fL->beam);
  for (i = 0; i < fL->n_archivos; i++) {
    for (j = 0; j < fL->fg_image->npixels; j++) {
      free (fL->atten[i][j + 1]);
    }
    free (fL->atten[i]);
  }
  free (fL->atten);

  if (fL->malla != NULL) {
    eliminarMallaVoronoi (fL->malla);
    fL->malla = NULL;
  }
  if (fL->imagen != NULL) {
    for (i = 0; i < fL->fg_image->size[0]; i++) {
      free (fL->imagen[i]);
    }
    free (fL->imagen);
    fL->imagen = NULL;
  }

  delete_map (fL->cmb_image);
  delete_map (fL->fg_image);

  for (i = 0; i < fL->fg_image->size[0]; i++) {
    free (fL->fourierI[i]);
    free (fL->mask[i]);
    free (fL->PixIntegrados_I[i]);
    free (fL->PixIntegrados_x[i]);
  }
  free (fL->fourierI);
  free (fL->mask);
  free (fL->PixIntegrados_I);
  free (fL->PixIntegrados_x);

  free (fL);
  fL = NULL;
}

funcL *copiarFuncL (funcL *fL_ori) {
  funcL *fL_copia;

  if (fL_ori == NULL) {
    return NULL;
  }

  fL_copia = newFuncL (fL_ori->nombreFits, fL_ori->infile, fL_ori->n_archivos, 
		       fL_ori->n_pols, fL_ori->entropia, -1, 
		       fL_ori->expanded_lx, fL_ori->expanded_ly, fL_ori->beam.fwhm / RPARCM);
  printf ("CoPiAndO con parametros %s, %s, %d, %d, %d, %g\n",
	  fL_ori->nombreFits, fL_ori->infile[0], fL_ori->n_archivos, 
	  fL_ori->n_pols, fL_ori->entropia, fL_ori->difmapNoise);
  if (fL_ori->malla != NULL) {
    NodoLista *n = fL_ori->malla->poligonos->primero;
    PoligonoVoronoi *pol;
    int i;

    printf ("Copiando Malla\n");
    fL_copia->malla = newMallaVoronoi();
    for (i = 0; i < fL_ori->malla->nPols - 3; i++) {
      pol = (PoligonoVoronoi *) n->info;

      insertarSitio (fL_copia->malla, &(pol->x), &(pol->y), pol->valor);

      n = n->sgte;
    }
  }
  printf ("copia : %d\n", fL_copia->n_pols);
  printf ("copia fL : %g\n", fL_copia);
  return fL_copia;
}

void normalizarVisibilidades (struct uvf_header **header,
			      struct uvf_sample **samples, int narch) {
  int iff, samp, arch;
  double visMax = -1e300, mod;

  for (arch = 0; arch < narch; arch++) {
    for (iff = 0; iff < header[arch]->nif; iff++) {
      for (samp = 0; samp < header[arch]->nsamp; samp++) {
	if ((mod = sqrt (SQR(samples[arch][samp].rdata[3 * iff]) +
			SQR(samples[arch][samp].rdata[3 * iff + 1]))) > visMax) {
	  visMax = mod;
	}
      }
    }
  }

  for (arch = 0; arch < narch; arch++) {
    for (iff = 0; iff < header[arch]->nif; iff++) {
      for (samp = 0; samp < header[arch]->nsamp; samp++) {
	//	printf ("Norm: %g + i %g\n", samples[arch][samp].rdata[3 * iff],
	//samples[arch][samp].rdata[3 * iff + 1]);
	samples[arch][samp].rdata[3 * iff]     /= visMax;
	samples[arch][samp].rdata[3 * iff + 1] /= visMax;
	samples[arch][samp].rdata[3 * iff + 2] *= visMax;
	//printf ("   vs %g + i %g\n", samples[arch][samp].rdata[3 * iff],
	//samples[arch][samp].rdata[3 * iff + 1]);
      }
    }
  }
}

void guardarFits (double **im, int **mask,
		  char *nombre, char *nombreFits) {
  char s[50];
  int i, j;
  struct image *imagen;
  
  imagen = do_read (nombreFits);
  if (im != NULL) {
    sprintf (s, "!%sIm.fits", nombre);
    for (i = 0; i < imagen->size[0]; i++) {
      for (j = 0; j < imagen->size[1]; j++) {
	imagen->pixels[i + j * imagen->size[0]] =
	  im[i][j];
      }
    }
    do_write_fits (imagen, s);


  }
  if (mask != NULL) {
    sprintf (s, "!%sMask.fits", nombre);
    for (i = 0; i < imagen->size[0]; i++) {
      for (j = 0; j < imagen->size[1]; j++) {
	imagen->pixels[i + j * imagen->size[0]] =
	  mask[i][j];
      }
    }
    do_write_fits (imagen, s);
  }

  delete_map (imagen);
}

void expanded_mockcbi(funcL *fL, int arch) {
  int i, nx, ny, n1x, n1y;
  //struct image *fg  = do_read ("expanded_fg.fits");
  //struct image *cmb = do_read ("expanded_cmb.fits");
  struct image *fg  = do_read (fL->nombreFits);
  struct image *cmb = do_read (fL->nombreFits);
  double lx = fL->expanded_lx, ly = fL->expanded_ly;

  n1x = round (fabs (lx / (fg->cdelt[0] * 60)));
  n1y = round (fabs(ly / (fg->cdelt[1] * 60)));
  nx = round (pow (2, floor (log (n1x) / log (2)) + 1));
  ny = round (pow (2, floor (log (n1y) / log (2)) + 1));
  //printf ("dx = %g, dx = %g\n", fg->cdelt[0], fg->cdelt[0] * 60);
  //printf ("nx1 = %d, ny1 = %d\n", n1x, n1y);
  //printf ("nx = %d, ny = %d\n", nx, ny);
  //exit(0);

  for (i = 0; i < fg->npixels; i++) {
    fg->pixels[i]  = fL->fg_image->pixels[i];
    cmb->pixels[i] = fL->cmb_image->pixels[i];
  }

  resize_map (fg, nx, ny, fg->crpix[0], fg->crpix[1]); 
  resize_map (cmb, nx, ny, cmb->crpix[0], cmb->crpix[1]);

  //do_write_fits (fg, "!expanded_fg.fits");
  //do_write_fits (cmb, "!expanded_cmb.fits");

  mockcbi_sub (cmb, fg, 
	       fL->header_obs[arch], fL->samples_obs[arch], 
	       fL->samples_mod[arch], fL->beam);
  delete_map (fg);
  delete_map (cmb);
}

void mockCBI (funcL *fL) {
  int i, j, arch, iff, samp;
  double Re, Im, x, y, kx, u, v;
  Status status;

  if(!fL->cmb_image) {
    printf("Error con cmb_image\n");
    exit(0);
  }
  for (arch = 0; arch < fL->n_archivos; arch++) {
    copiar_uvf_samples(fL->samples_obs[arch], fL->samples_mod[arch], fL->header_obs[arch]->nsamp);
    status = do_clear("clear", fL->header_obs[arch], fL->samples_mod[arch]);   /* Visibilidades = 0 */
    if (status == FAILURE)
      exit(1);
    
    for (iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
      for (samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	visMod (fL, arch, iff, samp, &Re, &Im);
	
	fL->samples_mod[arch][samp].rdata[3 * iff] = Re; 
	fL->samples_mod[arch][samp].rdata[3 * iff + 1] = Im;
	fL->samples_mod[arch][samp].rdata[3 * iff + 2] = 
	  fL->samples_obs[arch][samp].rdata[3 * iff + 2];
      }
    }
  }
}

void visMod (funcL *fL, int arch, int iff, int samp, double *visR, double *visI) {
  int i, j;
  double x, y, u, v , kx;
  (*visR) = 0;
  (*visI) = 0;

  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      y = (j - fL->y0[arch]) * fL->dy[arch];
      x = (i - fL->x0[arch]) * fL->dx[arch];
      u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
      v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
      kx = (u * x + v * y);

      (*visR) += fL->fg_image->pixels[i + j * fL->fg_image->size[0]]
	* fabs(fL->dx[arch] * fL->dy[arch])
	* fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff] * cos (2 * PI * kx)
	* fL->fg_scale; // sqrt (1 - x*x - y*y);

      (*visI) += fL->fg_image->pixels[i + j * fL->fg_image->size[0]]
	* fabs(fL->dx[arch] * fL->dy[arch])
	* fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff] * sin (2 * PI * kx)
	* fL->fg_scale; // sqrt (1 - x*x - y*y);
    }
  }
}

/*
funcL *copiar (funcL *fL) {
  int i;
  funcL fLcopia;
  if (fL == NULL) {
    return NULL;
  }

  fLcopia = newFuncL (fL->nombreFits, fL->infile, fL->n_archivos, 
		      fL->n_pols, fL->entropia, fL->difmapNoise);

  if (fL->malla != NULL) {
    // Copiamos la malla.
    NodoLista *n;
    PoligonoVoronoi *pol;

    fLcopia->malla = newMallaVoronoi ();
    n = (NodoLista *) fL->malla->poligonos->primero;
    for (i = 0; i < fL->malla->nPols; i++) {
      if (n == NULL) {
	fprintf (stderr, "ERROR al copiar malla, n == NULL\n");
	exit (0);
      }
      pol = (PoligonoVoronoi *) n;
      insertarSitio (fLcopia->malla, pol->x, pol->y, pol->valor);
      n = n->sgte;
    }
  }

  return fLcopia;
}
*/
