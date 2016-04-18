#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "funcL.h"

#define ENTROPIA 0

funcL * newFuncL (char * nombreImagen, char **nombresVis, int nVis, int n) {
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
  fL->ref_freq = 30.0e9; // 31.5e9;
  //fL->difmapNoise = 19.9451;  // [mJy / beam] IC443
  fL->difmapNoise = 19.9451;  // [mJy / beam] punt
  //fL->difmapNoise = 11.0664;  // [mJy / beam] IC443
  //fL->difmapNoise = 9.29467;  // [mJy / beam] L1622 IRAS12
  //fL->difmapNoise = 16.8758;  // [mJy / beam] L1622 IRAS12
  //fL->difmapNoise = 7.41446;  // [mJy / beam] L1622_simu

  fL->chi2 = 0;
  fL->S = 0;
  fL->imagen = NULL;
  fL->cmb_image = do_read (nombreImagen);
  fL->fg_image = do_read (nombreImagen);
  fL->header_obs = NULL;
  fL->samples_obs = NULL;
  fL->samples_mod = NULL;
  init_beam (&(fL->beam));

  fL->beam.type = CBI;
  printf ("Usando beam CBI.\n");
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
  }

  xf = (PLANCK_H * fL->ref_freq) / (BOLTZ_K * TCMB);
  g = SQR(xf)*exp(xf)/SQR(exp(xf)-1.0);

  //pixel = fabs(fL->fg_image->cdelt[0] * RPDEG * fL->fg_image->cdelt[1] * RPDEG);
  /* pixel solid angle in ster */
  //L->fg_scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * pixel * SQR(fL->ref_freq);
  /* K to Jy pixel^-1 */
  fL->fg_scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * SQR(fL->ref_freq);
  /* K to Jy srad^-1 */

  bmaj = fL->fg_image->bmaj * RPDEG; // en radianes.
  bmin = fL->fg_image->bmin * RPDEG; // en radianes.
  fL->difmapNoise *= (4 * log(2) * 1e-3 / (PI * bmaj * bmin)); // [Jy / sr]
  printf ("difmap noise 1 = %g [Jy / sr] \n", fL->difmapNoise);
  fL->difmapNoise /= fL->fg_scale; // en Kelvin.
  printf ("difmap noise 2 = %g K\n", fL->difmapNoise);
  printf ("fg_scale = %g\n", fL->fg_scale);

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

  do_write_fits (fL->fg_image, "!aux.fits");
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
      fprintf (stderr, "ERROR en newFuncL, fL-?PixIntegrados[%d] = NULL.\n", i);
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
double L (funcL *fL, double *pars, int n) {
  int i, j, arch, samp, iff;
  double ret = 0, datamin = 1e300, datamax = -1e300;
  
  // comienza comprobacion
  /*
  double a = 1, b = 5, c = 9;
  for (i = 0; i < n; i++) {
    ret += SQR(a * pars [3 * i] + b * pars [3 * i + 1] + c * pars [3 * i + 2]);
  }
  return ret;
  */
  // finaliza comprobacion

  if (fL->malla != NULL) {
    eliminarMallaVoronoi(fL->malla);
    fL->malla = NULL;
  }
  
  fL->malla = newMallaVoronoi();
  for (i = 0; i < n; i++) {
    if (pars[3 * i] < 0) {
      pars[3 * i] = 0;
    }
    else if (pars[3 * i] > 1) {
      pars[3 * i] = 1;
    }
    if (pars[3 * i + 1] < 0) {
      pars[3 * i + 1] = 0;
    }
    else if (pars[3 * i + 1] > 1) {
      pars[3 * i + 1] = 1;
    }
    /*
    printf ("insertando sitio (%g, %g, %g)\n", 
	    pars[3 * i], pars[3 * i + 1], 
	    pars[3 * i + 2] * fL->difmapNoise);
    */
    insertarSitio (fL->malla, pars[3 * i], pars[3 * i + 1], 
		   pars[3 * i + 2] * fL->difmapNoise);
    if (datamin > pars[3 * i + 2] * fL->difmapNoise) {
      datamin = pars[3 * i + 2] * fL->difmapNoise;
    }
    if (datamax < pars[3 * i + 2] * fL->difmapNoise) {
      datamax = pars[3 * i + 2] * fL->difmapNoise;
    }
  }

  fL->fg_image->datamax = datamax;
  fL->fg_image->datamin = datamin;
  /*
  for (i = 0, nodo = fL->malla->poligonos->primero; 
       i < n; i++, nodo = nodo->sgte) {
    pol = (PoligonoVoronoi *) nodo->info;
    pol->valor = pars[n - 1 - i];
  }
  */
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

  for (arch = 0; arch < fL->n_archivos; arch++) {
    mockcbi_sub (fL->cmb_image, fL->fg_image, 
		 fL->header_obs[arch], fL->samples_obs[arch], 
		 fL->samples_mod[arch], fL->beam);
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

  
  if (ENTROPIA) {
    ret = ret / 2 - S (fL, pars, n);
  }
  else {
    ret = ret / 2;
  }
  //printf ("Chi2 = %g\n", ret);
  if (fL->imagen != NULL) {
    for (i = 0; i < fL->fg_image->size[0]; i++) {
      free (fL->imagen[i]);
    }
    free (fL->imagen);
    fL->imagen = NULL;
  }
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

void dL (funcL *fL, double *pars, double *grad, int n, int aprox){
  int i, j, t, arch;
  double datamin = 1e300, datamax = -1e300;
  
  // comienza comprobacion
  /*
  double a = 1, b = 5, c = 9;
  for (i = 0; i < n; i++) {
    grad [3 * i]     = 2 * a * (a * pars [3 * i] + b * pars [3 * i + 1] + 
			    c * pars [3 * i + 2]);
    grad [3 * i + 1] = 2 * b * (a * pars [3 * i] + b * pars [3 * i + 1] + 
				c * pars [3 * i + 2]);
    grad [3 * i + 2] = 2 * c * (a * pars [3 * i] + b * pars [3 * i + 1] + 
				c * pars [3 * i + 2]);
  }
  return;
  */
  // finaliza comprobacion

  if (fL->malla != NULL){
    eliminarMallaVoronoi(fL->malla);
    fL->malla = NULL;
  }
  
  fL->malla = newMallaVoronoi();

  for (i = 0; i < n; i++) {
    if (pars[3 * i] < 0) {
      pars[3 * i] = 0;
    }
    else if (pars[3 * i] > 1) {
      pars[3 * i] = 1;
    }

    if (pars[3 * i + 1] < 0) {
      pars[3 * i + 1] = 0;
    }
    else if (pars[3 * i + 1] > 1) {
      pars[3 * i + 1] = 1;
    }
    
    insertarSitio (fL->malla, pars[3 * i], pars[3 * i + 1],
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
  if (ENTROPIA) {
    dS(fL, pars, grad, n);
  }
  fL->fg_image->datamax = datamax;
  fL->fg_image->datamin = datamin;
  
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

  for (arch = 0; arch < fL->n_archivos; arch++) {
    mockcbi_sub(fL->cmb_image, fL->fg_image, 
		fL->header_obs[arch], fL->samples_obs[arch],
		fL->samples_mod[arch], fL->beam);
  }
  
  printf ("Calculando dLdI\n");
  t = time (0);
  dLdI (fL, grad, n, aprox);
  printf ("tiempo dLdI = %d\n", (int) time (0) - t);
  printf ("Calculando dLdx\n");
  t = time(0);
  dLdx (fL, grad, n);
  printf ("tiempo dLdx = %d\n", (int) time (0) - t);
  printf ("Calculado dL\n");
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
    grad[3 * i + 2] = dS_cte;
    for (j = 1; j <= pars[3 * i + 2]; j++) {
      grad[3 * i + 2] += 1.0 / j;
    }
  }
}

/* Calculo del Gradiente con respecto a las posiciones de los
 * poligonos. 
 */

void dLdx (funcL *fL, double *grad, int n){
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
	  dVdx (fL, arch, iff, samp, pol, &dVdx_R, &dVdx_I, &dVdy_R, &dVdy_I);
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

// Calculo del Gradiente con respecto a I.

void dLdI (funcL *fL, double *grad, int n, int aprox) {
  int i, j, arch, idPol;

  printf ("Calculando Fourier\n");
  if (aprox) {
    calcularFourierIAprox (fL);
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
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      if (fL->mask[i][j] < 0) {
	fprintf (stderr, "ERROR en dLdI: fL->mask[%d][%d] = %d\n", 
		 i, j, fL->mask[i][j]);
	exit (1);
      }
      idPol = fL->mask[i][j] - 3;
      grad[3 * idPol + 2] += fL->fourierI[i][j] * fL->difmapNoise * fL->fg_scale;
    }
  }
  
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
  
void calcularSin (funcL *fL) {
  int arch, samp, iff;
  double dx, dy;//, x, y;
  double u, v, sinu, sinv;

  //dx = fL->fg_image->cdelt[0] * 60;   // en arcmin
  //dy = fL->fg_image->cdelt[1] * 60;   // en arcmin
  //dx = fL->fg_image->cdelt[0] * RPDEG;   // en radianes
  //dy = fL->fg_image->cdelt[1] * RPDEG;   // en radianes
  dx = fabs (fL->fg_image->cdelt[0] * RPDEG);   // en radianes GUILLE
  dy = fabs (fL->fg_image->cdelt[1] * RPDEG);   // en radianes GUILLE

  printf ("dx = %g, dy = %g\n", dx, dy);

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
	  sinu = sin (PI * u * dx) / (PI * u);
	}
	else {
	  sinu = dx;
	  //sinu = dx / RPDEG;
	  printf ("sinu 0, samp = %d, iff = %d, u = %g, iffreq = %g\n", 
		  samp, iff, fL->samples_obs[arch][samp].u, 
		  fL->header_obs[arch]->iffreq[iff]);
	}
	if (v != 0) {
	  sinv = sin (PI * v * dy) / (PI * v);
	}
	else {
	  sinv = dy;
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
    //dx = fL->fg_image->cdelt[0] * RPDEG;   // en radianes
    //dy = fL->fg_image->cdelt[1] * RPDEG;   // en radianes
    dx = fabs (fL->fg_image->cdelt[0] * RPDEG);   // en radianes GUILLE
    dy = fabs (fL->fg_image->cdelt[1] * RPDEG);   // en radianes GUILLE
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
	   double *dVdx_R, double *dVdx_I, double *dVdy_R, double *dVdy_I) {
  AristaVoronoi *a;
  double integralR = 0, integralI = 0;
  
  (*dVdx_R) = 0;
  (*dVdx_I) = 0;
  (*dVdy_R) = 0;
  (*dVdy_I) = 0;
  //a = (pol->a->poliDer == pol)? pol->a->cwPred : pol->a->cwSucc;
  a = pol->a;
  do {
    integraldVdx (fL, arch, iff, samp, pol, a, &integralR, &integralI);
    //integralR *= fL->fg_scale;
    //integralI *= fL->fg_scale;
    if (a->poliDer == pol) {
      (*dVdx_R) += (pol->valor - a->poliIzq->valor) * a->cosDer * integralR;
      (*dVdy_R) += (pol->valor - a->poliIzq->valor) * a->sinDer * integralR;
      (*dVdx_I) += (pol->valor - a->poliIzq->valor) * a->cosDer * integralI;
      (*dVdx_I) += (pol->valor - a->poliIzq->valor) * a->sinDer * integralI;
      a = a->cwPred;
    }
    else if (a->poliIzq == pol) {
      (*dVdx_R) += (pol->valor - a->poliDer->valor) * a->cosIzq * integralR;
      (*dVdy_R) += (pol->valor - a->poliDer->valor) * a->sinIzq * integralR;
      (*dVdx_I) += (pol->valor - a->poliDer->valor) * a->cosIzq * integralI;
      (*dVdx_I) += (pol->valor - a->poliDer->valor) * a->sinIzq * integralI;
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
  (*dVdx_R) *= fL->fg_scale * 0.5;
  (*dVdx_I) *= fL->fg_scale * 0.5;
  (*dVdy_R) *= fL->fg_scale * 0.5;
  (*dVdy_I) *= fL->fg_scale * 0.5;
}

void integraldVdx (funcL *fL, int arch, int iff, int samp, 
		   PoligonoVoronoi *pol, AristaVoronoi *a,
		   double *integralR, double *integralI) {
  int i, j, nx, ny, tipoLinea, lineaAnterior = 0;
  PuntoVoronoi *p1, *p2, *pIni, *pFin, *paux;
  double m, u, v, I1R, I1I, c1, c2;
  double x1, y1, x2, y2, s0, sin_alfa, cos_alfa, t, dt, sin_c1dt;
  
  p1 = newPuntoVoronoi (0, 0);
  p2 = newPuntoVoronoi (0, 0);
  nx = fL->fg_image->size[0];
  ny = fL->fg_image->size[1];
  (*integralR) = 0;
  (*integralI) = 0;
  u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
  v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
  
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
      return;
    }
  }
  else {
    if (!interseccionCuadrado (a->ptoFin, a->ptoIni, p2, p1, tipoLinea)) {
      free (p1);
      free (p2);
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
  //x1 = (pIni->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
  //y1 = (pIni->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
  x1 = (p1->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
  y1 = (p1->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
  s0 = y1 * sin_alfa + x1 * cos_alfa;
  c1 = PI * (v * cos_alfa - u * sin_alfa);
  c2 = PI * s0 * (u * cos_alfa + v * sin_alfa);

  pFin = newPuntoVoronoi (p1->x, p1->y);
  
  for (pIni = newPuntoVoronoi (p1->x, p1->y);
       (pFin->x != p2->x) && (pFin->y != p2->y);) {
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
	SQR(1.0/(nx - 1.0)) + SQR(1.0/(ny - 1.0))) {
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
      exit (1);    
      }
    
    x1 = (pIni->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
    y1 = (pIni->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
    x2 = (pFin->x * (nx - 1.0) - fL->x0[arch]) * fL->dx[arch];
    y2 = (pFin->y * (ny - 1.0) - fL->y0[arch]) * fL->dy[arch];
    
    t = 0.5 * (cos_alfa * (y1 + y2) - sin_alfa * (x1 + x2));
    dt = sin_alfa * (x1 - x2) + cos_alfa * (y2 - y1);
    sin_c1dt = sin (c1 * dt);
    // Trans Fourier
    if (c1 != 0) {
      I1R = cos (2 * (c2 + c1 * t)) * sin_c1dt / c1;
      I1I = sin (2 * (c2 + c1 * t)) * sin_c1dt / c1;
    } else {
      I1R = cos (2 * (c2 + c1 * t));
      I1I = sin (2 * (c2 + c1 * t));
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
    fL->PixIntegrados_x[i][j]++;
    
    pIni->x = pFin->x;
    pIni->y = pFin->y;
  }
  free (pIni);
  free (pFin);
  
  free (p1);
  free (p2);
}

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
  printf ("lineaAnterior = %d, (i, j) = (%d, %d)\n", lineaAnterior, i, j);
  printf ("p1 = (%g, %g), pFin = (%g, %g)\n", p1->x, p1->y, pFin->x, pFin->y);
  printf ("p1 = (%g, %g), pFin = (%g, %g) (en pixeles)\n", 
	  p1->x * (nx - 1), p1->y * (ny - 1), 
	  pFin->x * (nx - 1), pFin->y * (ny - 1));
  printf ("p1 - (i, j) = (%g, %g)\n", 
	  p1->x * (nx - 1), p1->y * (ny - 1) - j - 0.5);  
  */
  // linea inferior.
  if (lineaAnterior != BORDE_SUP) {
    y = (j - 0.5) / (ny - 1.0);
    x = (y - p1->y) * dx1 / dy1 + p1->x;
    //printf ("inferior y = %g, %g en pixeles, ppunto = %g\n", y, y * (ny - 1),
    //(x - p1->x) * (pFin->x - p1->x) + 
    //(y - p1->y) * (pFin->y - p1->y));
    if ((x - p1->x) * (pFin->x - p1->x) + 
	(y - p1->y) * (pFin->y - p1->y) > 0) {
      //d1 = pow (p1->x - x, 2) + pow (p1->y - y, 2);
      d1 = SQR (p1->x - x) + SQR (p1->y - y);
      //printf ("OK1 d1 = %g, d = %g\n", d1, d);
      d = d1;
      pSgte->x = x;
      pSgte->y = y;
      ret = BORDE_INF;
      //printf ("OK, (%g, %g) (en pixeles)\n", x * (nx - 1), y * (ny - 1));
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
	(y - p1->y) * (pFin->y - p1->y) > 0) {
      //d1 = pow (p1->x - x, 2) + pow (p1->y - y, 2); 
      d1 = SQR (p1->x - x) + SQR (p1->y - y); 
      //printf ("OK1 d1 = %g, d = %g\n", d1, d);
      if (d1 < d) {
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
	(y - p1->y) * (pFin->y - p1->y) > 0) {
      //d1 = pow (p1->x - x, 2) + pow (p1->y - y, 2); 
      d1 = SQR (p1->x - x) + SQR (p1->y - y); 
      //printf ("OK1 d1 = %g, d = %g\n", d1, d);
      if (d1 < d) {
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
    //    (x - p1->x) * (pFin->x - p1->x) + 
    //    (y - p1->y) * (pFin->y - p1->y));
    //printf ("(x, y) = (%g, %g)\n", x * (nx - 1), y * (ny - 1));
    if ((x - p1->x) * (pFin->x - p1->x) + 
	(y - p1->y) * (pFin->y - p1->y) > 0) {
      //d1 = pow (p1->x - x, 2) + pow (p1->y - y, 2); 
      d1 = SQR (p1->x - x) + SQR (p1->y - y); 
      //printf ("OK1 d1 = %g, d = %g\n", d1, d);
      if (d1 < d) {
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
  if (d >= d2) {
    pSgte->x = pFin->x;
    pSgte->y = pFin->y;
    //printf ("OK Fin, d = %g, d2 = %g\n", d, d2);
  }
  //printf ("ultimo OK\n\n");
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
	for (iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
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
	    //fL->fourierI[i][j] += fL->sin[arch][samp][iff] 29 Oct 2005
	    //* fL->samples_mod[arch][samp].rdata[3 * iff + 2] * (Re + Im);
	    fL->fourierI[i][j] += fL->sin[arch][samp][iff] 
	      * fL->samples_mod[arch][samp].rdata[3 * iff + 2] * (Re + Im)
	      * fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff];
	  }
	  //fL->fourierI[i][j] *= fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff]; 29 Oct 2005
	}
      }
    }
  }
}

void calcularFourierIAprox (funcL *fL) {
  int i, j, arch, iff, samp;
  double Re, Im, x, y, kx, u, v, dx, dy;
  
  //dx = fL->fg_image->cdelt[0] * RPDEG;   // en radianes
  //dy = fL->fg_image->cdelt[1] * RPDEG;   // en radianes
  dx = fabs (fL->fg_image->cdelt[0] * RPDEG);   // en radianes GUILLE
  dy = fabs (fL->fg_image->cdelt[1] * RPDEG);   // en radianes GUILLE
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
	  }
	  fL->fourierI[i][j] *= fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff];
	}
      }
      fL->fourierI [i][j] *= dx * dy;
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
