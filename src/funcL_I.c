#include <stdio.h>
#include <math.h>
#include "funcL.h"

funcL * newFuncL (char * nombreImagen, char **nombresVis, int nVis, int n){
  int i, j;
  double xf, g;//, pixel;
  Status status;
  funcL *L = malloc (sizeof (funcL));

  L->malla = newMallaVoronoi();
  for (i = 0; i < n; i++) {
    insertarSitio (L->malla, (double) random() / RAND_MAX,
		   (double) random() / RAND_MAX, 0.0);
  }

  L->n_pols = n;
  L->ref_freq = 30.0e9; // 31.5e9;
  L->chi2 = 0;
  L->S = 0;
  L->imagen = NULL;
  L->cmb_image = do_read(nombreImagen);
  L->fg_image = do_read(nombreImagen);
  L->header_obs = NULL;
  L->samples_obs = NULL;
  L->samples_mod = NULL;
  init_beam(&(L->beam));

  L->beam.type = CBI;
  printf("Usando beam CBI.\n");
  // beam.cutoff *= 0.1;

  if(!L->cmb_image)
  {
    fprintf(stderr, "ERROR al intentar abrir %s.", nombreImagen);
    exit(0);
  }

  if (!L->cmb_image->bunit || !(strcmp("K", L->cmb_image->bunit) == 0))
  {
    fprintf(stderr, "ERROR: Unidades de %s  no estan en K. bunit = %s\n", 
	    nombreImagen, L->cmb_image->bunit);
    //exit (1);
  }
  if(!L->fg_image)
  {
    fprintf(stderr, "ERROR al intentar abrir %s.", nombreImagen);
    exit(1);
  }
  xf = (PLANCK_H * L->ref_freq) / (BOLTZ_K * TCMB);
  g = SQR(xf)*exp(xf)/SQR(exp(xf)-1.0);
  //pixel = fabs(L->fg_image->cdelt[0] * RPDEG * L->fg_image->cdelt[1] * RPDEG);
  /* pixel solid angle in ster */
  //L->fg_scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * pixel * SQR(L->ref_freq);
  /* K to Jy pixel^-1 */
  L->fg_scale = (2e26*BOLTZ_K/SQR(LIGHTSPEED)) * g * SQR(L->ref_freq);
  /* K to Jy srad^-1 */

  L->n_archivos = nVis;
  L->header_obs = (struct uvf_header**) malloc(L->n_archivos * 
					    sizeof(struct uvf_header *));
  L->samples_obs = (struct uvf_sample **) malloc(L->n_archivos * 
					     sizeof(struct uvf_sample *));
  L->samples_mod = (struct uvf_sample **) malloc(L->n_archivos * 
					     sizeof(struct uvf_sample *));
  //printf("ok2\n");
  L->infile = malloc(L->n_archivos * sizeof(char *));
  for(i = 0; i < L->n_archivos; i++)
  // abrimos los archivos de visibilidades
    {
      status = do_read_uvdata(nombresVis[i], &(L->header_obs[i]), 
			      &(L->samples_obs[i]));
      //multiplicar_vis(L->header_obs[i], L->samples_obs[i], norm_factor);
      if(status == SUCCESS)
	{
	  L->infile[i] = newstring(nombresVis[i]);
	}
      else
	{
	  printf("Error con el archivo uvf\n");
	  exit(1);
	}
    }
   for (i = 0; i < L->n_archivos; i++) {
      L->samples_mod[i] = (struct uvf_sample*) malloc((L->header_obs[i]->nsamp) * 
						    sizeof(struct uvf_sample));
   }
  // ya tenemos abiertos los archivos .sub
  
  L->atten = atenuacion(L->fg_image, L->header_obs, L->n_archivos, L->beam);

  L->mask = (int **) malloc(L->fg_image->size[0] * sizeof(int *));
  L->PixIntegrados = (int **) malloc(L->fg_image->size[0] * sizeof(int *));
  for (i = 0; i < L->fg_image->size[0]; i++) {
    L->mask[i] = (int *) malloc (L->fg_image->size[1]  * sizeof(int));
    L->PixIntegrados[i] = (int *) malloc (L->fg_image->size[1]  * sizeof(int));
    for (j = 0; j < L->fg_image->size[1]; j++) {
      L->PixIntegrados[i][j] = 0;
    }
  }

  printf ("Calculando sincos\n");
  calcularSin (L);
  calculardxdy (L);
  printf ("L listo\n");

  L->nPixIntegrados = 0;
  return L;
}

/* Calcula la funcion L.
 * En pars los iesimos parametros son la posicion x, los i+1 la posicion y
 * y los i+2 las intensidades de los n centros de Voronoi.
 */
double L (funcL *fL, double *pars, int n){
  int i, j, k;
  double ret = 0, S = 0;
  //double *im;
  NodoLista *nodo;
  PoligonoVoronoi *pol;
  /*
  fL->malla = newMallaVoronoi();
  for (i = 0; i < n; i++) {
    insertarSitio (fL->malla, pars[3 * i], pars[3 * i + 1], pars[3 * i + 2]);
  }
  */
  for (i = 0, nodo = fL->malla->poligonos->primero; 
       i < n; i++, nodo = nodo->sgte) {
    pol = (PoligonoVoronoi *) nodo->info;
    pol->valor = pars[n - 1 - i];
    //printf ("pol %d = %f vs pars[%d] = %g\n", pol->id, (float) pol->valor, 
    //    n - 1 - i, pars[n - 1 - i]);
  }
  //im = malloc (fL->fg_image->npixels * sizeof (double));
  fL->imagen = toImage(fL->malla, fL->fg_image->size[0], fL->fg_image->size[1], 
		       fL->mask);
  
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      fL->fg_image->pixels[i + j * fL->fg_image->size[0]] = fL->imagen[i][j];
    }
  }
  for (k = 0; k < fL->n_archivos; k++) {
    mockcbi_sub(fL->cmb_image, fL->fg_image, 
		fL->header_obs[k], fL->samples_obs[k], fL->samples_mod[k], 
		fL->beam);
    for(j = 0; j < fL->header_obs[k]->nchan; j++) {
      for(i = 0; i < fL->header_obs[k]->nsamp; i++) {
	ret += ((SQR(fL->samples_mod[k][i].rdata[j * 3] 
		     - fL->samples_obs[k][i].rdata[j * 3])
		 + SQR(fL->samples_mod[k][i].rdata[j * 3 + 1] 
		       - fL->samples_obs[k][i].rdata[j * 3 + 1]))
		* fL->samples_mod[k][i].rdata[j * 3 + 2]);
      }
    }
  }
  ret = ret / 2 - S;
  //printf ("Chi2 = %g\n", ret);
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    free (fL->imagen[i]);
  }
  free (fL->imagen);
  //free (im);
  //eliminarMallaVoronoi(fL->malla);
  /*  for (k = 0; k < fL->n_archivos; k++) {
    free(fL->samples_mod[k]);
    }*/
  return ret;
}

void dL (funcL *fL, double *pars, double *grad, int n){
  int i, j, arch, chan, samp;
  double dVdI_R, dVdI_I, Re, Im;
  PoligonoVoronoi *pol;
  NodoLista *nodo;

  /*
  L->malla = newMallaVoronoi();
  for (i = 0; i < n; i++) {
    insertarSitio (L->malla, pars[3 * i], pars[3 * i + 1], pars[3 * i + 2]);
  }
  */

  fL->nPixIntegrados = 0;
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      fL->PixIntegrados[i][j] = 0;
    }
  }

  for (i = 0, nodo = fL->malla->poligonos->primero; 
       i < n; i++, nodo = nodo->sgte) {
    pol = (PoligonoVoronoi *) nodo->info;
    /*
      if ((pol->x != pars[3 * i]) || (pol->y != pars[3 * i + 1])) {
      fprintf (stderr, "ERROR: dL lista no coincide con arreglo. :P\n");
      exit (1);
      }
    */
    if (pol->id - 3 != n - 1 - i) {
      fprintf (stderr, "ERROR dL: pol->id = %d != %d\n", pol->id, n - 1 - i);
      exit (1);
    }
    pol->valor = pars[n - 1 - i];
  }

  fL->imagen = toImage(fL->malla, fL->fg_image->size[0], fL->fg_image->size[1], 
		       fL->mask);
  
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[0]; j++) {
      //printf ("fL->imagen[%d][%d] = %g\n", i, j, fL->imagen[i][j]);
      fL->fg_image->pixels[i + j * fL->fg_image->size[0]] = fL->imagen[i][j];
    }
  }
  
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    free (fL->imagen[i]);
  }
  free (fL->imagen);
  
  for (arch = 0; arch < fL->n_archivos; arch++) {
    mockcbi_sub(fL->cmb_image, fL->fg_image, 
		fL->header_obs[arch], fL->samples_obs[arch],
		fL->samples_mod[arch], fL->beam);
    
    //for (i = 0, nodo = fL->malla->poligonos->primero; 
    //i < 156; i++, nodo= nodo->sgte);
    for (i = 0, nodo = fL->malla->poligonos->primero; 
     i < n; i++, nodo = nodo->sgte) {
      pol = (PoligonoVoronoi *) nodo->info;
      //printf ("pol%d = (%g, %g) = %g", pol->id, pol->x, pol->y, pol->valor);
      //grad[i] = 0;
      grad[n - 1 - i] = 0;
      Re = 0;
      Im = 0;
      for (chan = 0; chan < fL->header_obs[arch]->nchan; chan++) {
	for (samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	  dVdI_R = 0;
	  dVdI_I = 0;
	  dVdI (fL, arch, chan, samp, pol, &dVdI_R, &dVdI_I);
	  Re = (fL->samples_mod[arch][samp].rdata[3 * chan] 
		- fL->samples_obs[arch][samp].rdata[3 * chan]) * dVdI_R;
	  Im = (fL->samples_mod[arch][samp].rdata[3 * chan + 1] 
		- fL->samples_obs[arch][samp].rdata[3 * chan + 1]) * dVdI_I;

	  //grad[i] += fL->samples_mod[arch][samp].rdata[3 * chan + 2] * (Re + Im);
	  grad[n - 1 - i] += 
	    fL->samples_mod[arch][samp].rdata[3 * chan + 2] * (Re + Im);
	  if (isnan(grad[n - 1 - i])) {
	    fprintf(stderr, "ERROR dL: grad[%d] = %g\n", i, grad[n - 1 - i]);
	    fprintf(stderr, "          Re = %g, Im = %g\n", Re, Im);
	    fprintf(stderr, "          dVR = %g, dVI = %g\n", dVdI_R, dVdI_I);
	    exit (1);
	  }
	}
      }
      //grad[i] *= fL->fg_scale;
      grad[n - 1 - i] *= fL->fg_scale ;
      //printf("   grad[%d] = %g\n", n - 1 - i, grad[n - 1 - i]);
    }
    if (fL->nPixIntegrados != 
	fL->fg_image->npixels * fL->header_obs[0]->nchan * 
	fL->header_obs[0]->nsamp) {
      struct image *pix_image;
      fprintf(stderr, "ERROR dL: Distinto numero de pixeles integrados\n");
      fprintf(stderr, "          nPixIntegrados = %d, npixels = %d\n", 
	      fL->nPixIntegrados, fL->fg_image->npixels);
      pix_image = do_read ("mask_0.fits");
      for (i = 0; i < pix_image->size[0]; i++) {
	for (j = 0; j < pix_image->size[1]; j++) {
	  pix_image->pixels[i + j * pix_image->size[0]] = 
	    fL->PixIntegrados[i][j];
	}
      }
      do_write_fits (pix_image, "!pixInt_ERROR.fits");
      
      exit (1);
    }
  }
}
  
void calcularSin(funcL *L){
  int arch, samp, chan;
  double dx, dy;//, x, y;
  double u, v, sinu, sinv;

  //dx = L->fg_image->cdelt[0] * 60;   // en arcmin
  //dy = L->fg_image->cdelt[1] * 60;   // en arcmin
  dx = L->fg_image->cdelt[0] * RPDEG;   // en radianes
  dy = L->fg_image->cdelt[1] * RPDEG;   // en radianes

  L->sin = (double ***) malloc (L->n_archivos * sizeof (double **));
  for (arch = 0; arch < L->n_archivos; arch++) {
    L->sin[arch] = (double **) malloc(L->header_obs[arch]->nsamp 
				      * sizeof(double*));
    
    for(samp = 0; samp < L->header_obs[arch]->nsamp; samp++) {
      L->sin[arch][samp] = (double *) malloc(L->header_obs[arch]->nchan 
					     * sizeof(double));
      for (chan = 0; chan < L->header_obs[arch]->nchan; chan++) {
	u = L->samples_obs[arch][samp].u * L->header_obs[arch]->iffreq[chan];
	v = L->samples_obs[arch][samp].v * L->header_obs[arch]->iffreq[chan];
	if (u != 0) {
	  sinu = sin (PI * u * dx) / (PI * u);
	}
	else {
	  sinu = dx / RPDEG;
	  printf ("sinu 0\n");
	}
	if (v != 0) {
	  sinv = sin (PI * v * dy) / (PI * v);
	}
	else {
	  sinv = dy / RPDEG;
	  printf ("sinv 0\n");
	}
	/*
	L->sin[arch][samp][chan] = (sin (PI * u * dx) * sin (PI * v * dy) 
				    / (PI * PI * u * v));
	*/
	L->sin[arch][samp][chan] = sinu * sinv;
	if (isnan(L->sin[arch][samp][chan])) { 
	  fprintf(stderr, "ERROR calcularSin: sin[%d][%d][%d] = %g\n", 
		  arch, samp, chan, L->sin[arch][samp][chan]);
	  exit (1);
	}
      }
    }
  }
}

void calculardxdy (funcL *L){
  double lobs, mobs, obsra, obsdec, raimage, decimage, ra_offset = 0, dx, dy;
  int i;

  raimage = L->fg_image->crval[0] * RPDEG;
  decimage = L->fg_image->crval[1] * RPDEG;

  L->x0 = (int *) malloc (L->n_archivos * sizeof (int));
  L->y0 = (int *) malloc (L->n_archivos * sizeof (int));
  L->dx = (double *) malloc (L->n_archivos * sizeof (double));
  L->dy = (double *) malloc (L->n_archivos * sizeof (double));

  for (i = 0; i < L->n_archivos; i++) {
    obsra = ra_offset + L->header_obs[i]->obsra * RPDEG;
    obsdec = L->header_obs[i]->obsdec * RPDEG;  
    direccos(obsra, obsdec, raimage, decimage, &lobs, &mobs);
    dx = L->fg_image->cdelt[0] * RPDEG;   // en radianes
    dy = L->fg_image->cdelt[1] * RPDEG;   // en radianes
    L->x0[i] = (L->fg_image->crpix[0]-1.0) + lobs/dx;
    L->y0[i] = (L->fg_image->crpix[1]-1.0) + mobs/dy;
    //L->dx[i] = L->fg_image->cdelt[0] * 60;   // en arcmin
    //L->dy[i] = L->fg_image->cdelt[1] * 60;   // en arcmin
    L->dx[i] = dx;   // en radianes
    L->dy[i] = dy;   // en radianes
  }
}

void dVdI (funcL *fL, int arch, int chan, int samp, 
	   PoligonoVoronoi *pol, double *dVdI_R, double *dVdI_I) {
  int i, j, i_ini, j_ini, i1, i2, i11 = -1, i21 = -1;

  i_ini = round(pol->x * (fL->fg_image->size[0] - 1));
  j_ini = round(pol->y * (fL->fg_image->size[1] - 1));
  if (pol->id != fL->mask[i_ini][j_ini]) {
    for (i = i_ini - 1;i <= i_ini + 1; i++) {
      for (j = j_ini - 1; j <= j_ini + 1; j++) {
	if (i >= 0 && i < fL->fg_image->size[0] && 
	    j >= 0 && j < fL->fg_image->size[1] &&
	    pol->id == fL->mask[i][j]) {
	  i_ini = i;
	  j_ini = j;
	  i = i_ini + 2;
	  j = j_ini + 2;
	}
      }
    }   
    if (pol->id != fL->mask[i_ini][j_ini]) {
      struct image *mask_image;
      mask_image = do_read ("mask_0.fits");
      fprintf (stderr, "ERROR dVdI: id = %d, mask[%d][%d] = %d\n", 
	       pol->id, i_ini, j_ini, fL->mask[i_ini][j_ini]);
      for (i = 0; i < mask_image->size[0]; i++) {
	for (j = 0; j < mask_image->size[1]; j++) {
	  mask_image->pixels[i + j * mask_image->size[0]] = fL->mask[i][j];
	}
      }
      do_write_fits(mask_image, "!mask_ERROR.fits");
      
      exit (1);
    }
  }
  (*dVdI_R) = 0;
  (*dVdI_I) = 0;
  i = i_ini;
  i2 = i + 1;

  for (j = j_ini; i2 + 1 >= i && j < fL->fg_image->size[1]; j++) {
    integraldVdIH (fL, arch, chan, samp, pol, dVdI_R, dVdI_I, i, j, &i1, &i2);
    if (i11 == -1) {
      i11 = i1;
      i21 = i2;
    }
    
    if (isnan(*dVdI_R) || isnan(*dVdI_I)) {
      fprintf (stderr, "ERROR dVdI: j = %d, dVR = %g, dVI = %g\n",
	      j, *dVdI_R, *dVdI_I);
      exit (1);
    }
    
    if (j + 1 < fL->fg_image->size[1]) {
      for (i = (i1 <= 0)? 0 :i1 - 1; 
	   i < fL->fg_image->size[0] && pol->id != fL->mask[i][j+1] && i <= i2 + 1;
	   i++);
    }
  }

  j = j_ini;
  if (j - 1 >= 0) {
    for (i = (i11 <= 0)? 0 :i11 - 1; 
	 i < fL->fg_image->size[0] && pol->id != fL->mask[i][j-1] && i <= i21 + 1;
	 i++);
  }
  i2 = i21;
  for (j = j_ini - 1; i2 + 1 >= i && j >= 0; j--) {
    integraldVdIH (fL, arch, chan, samp, pol, dVdI_R, dVdI_I, i, j, &i1, &i2);
    
    if (isnan(*dVdI_R) || isnan(*dVdI_I)) {
      fprintf(stderr, "ERROR dVdI: j = %d, dVR = %g, dVI = %g\n", 
	      j, *dVdI_R, *dVdI_I);
      exit (1);
    }
    
    if (j - 1 >= 0) {
      for (i = (i1 <= 0)? 0 :i1 - 1; 
	   i < fL->fg_image->size[0] && pol->id != fL->mask[i][j-1] && i <= i2 + 1;
	   i++);
    }
  }
}

/* Integra horizontalmenete a la altura j, dejando en i1 e i2 el primer 
 * y ultimo punto sumado respectivamente.
 */
void integraldVdIH (funcL *fL, int arch, int chan, int samp, 
		    PoligonoVoronoi *pol, double *dVdI_R, double *dVdI_I,
		    int i_ini, int j_ini, int *i1, int *i2) {
  
  double x, y, factor, kx, u, v;
  int i, j;
  
  j = j_ini;
  
  if ((j < 0) || (j_ini >= fL->fg_image->size[1])
      || (i_ini < 0) || (i_ini >= fL->fg_image->size[0])) {
    return;
  }

  if (pol->id != fL->mask[i_ini][j_ini]){
    return;
  }
  
  u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[chan];
  v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[chan];
  y = (j - fL->y0[arch]) * fL->dy[arch];
  
  i = i_ini;
  if (i >= fL->fg_image->size[0]) {
    i = fL->fg_image->size[0] - 1;
  }
  for (; i >= 0 && pol->id == fL->mask[i][j]; i--) { 
    x = (i - fL->x0[arch]) * fL->dx[arch];
    factor = (fL->atten[arch][i + j * fL->fg_image->size[0] + 1][chan]
	      * fL->sin[arch][samp][chan]);
    kx = (u * x + v * y);
    
    (*dVdI_R) += factor * cos (2 * PI * kx);
    (*dVdI_I) += factor * sin (2 * PI * kx);
    fL->nPixIntegrados++;
    fL->PixIntegrados[i][j]++;
    
    if (isnan(*dVdI_R) || isnan(*dVdI_I)) { 
      fprintf(stderr, "ERROR integraldVdIH: i = %d, dVR = %g, dVI = %g\n", 
	      j, *dVdI_R, *dVdI_I);
      fprintf(stderr, "                     factor = %g\n", factor);
      fprintf(stderr, "                     atten = %g, sin = %g \n", 
	      fL->atten[arch][i + j * fL->fg_image->size[0] + 1][chan],
	      fL->sin[arch][samp][chan]);
      exit (1);
    }
    
  }
  (*i1) = i + 1;
  
  i = i_ini + 1;
  if (i < 0) {
    i = 0;
  }
  for (; i < fL->fg_image->size[0] && pol->id == fL->mask[i][j]; i++) { 
    x = (i - fL->x0[arch]) * fL->dx[arch];
    factor = (fL->atten[arch][i + j * fL->fg_image->size[0] + 1][chan]
	      * fL->sin[arch][samp][chan]);
    kx = (u * x + v * y);
    
    (*dVdI_R) += factor * cos (2 * PI * kx);
    (*dVdI_I) += factor * sin (2 * PI * kx);
    fL->nPixIntegrados++;
    fL->PixIntegrados[i][j]++;
  }
  (*i2) = i - 1;
}

void eliminarFuncL (funcL *L){
  int i, j;

  for (i = 0; i < L->n_archivos; i++) {
    free(L->samples_obs[i]);
    free(L->header_obs[i]);
    free(L->samples_mod[i]);
    free (L->infile[i]);
    L->infile[i] = NULL;
  }
  free (L->header_obs);
  free (L->samples_obs);
  free (L->samples_mod);
  free (L->infile);
  
  L->atten = atenuacion(L->fg_image, L->header_obs, L->n_archivos, L->beam);
  for (i = 0; i < L->n_archivos; i++) {
    for (j = 0; j < L->fg_image->npixels; j++) {
      free (L->atten[i][j + 1]);
    }
    free (L->atten[i]);
  }
  free (L->atten);

  delete_map (L->cmb_image);
  delete_map (L->fg_image);

  for (i = 0; i < L->fg_image->size[0]; i++) {
    free (L->mask[i]);
  }
  free (L->mask);
  
  free (L);
}
