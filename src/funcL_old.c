#include <stdio.h>
#include <math.h>
#include "funcL.h"

funcL * newFuncL (char * nombreImagen, char **nombresVis, int nVis){
  int i;
  Status status;

  funcL *L = malloc (sizeof (funcL));
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
    fprintf(stderr, "ERROR: Unidades de %s  no estan en K.\n", nombreImagen);
  }
  if(!L->fg_image)
  {
    fprintf(stderr, "ERROR al intentar abrir %s.", nombreImagen);
    exit(1);
  }
  L->n_archivos = nVis;
  L->header_obs = (struct uvf_header**) malloc(L->n_archivos * 
					    sizeof(struct uvf_header *));
  L->samples_obs = (struct uvf_sample **) malloc(L->n_archivos * 
					     sizeof(struct uvf_sample *));
  L->samples_mod = (struct uvf_sample **) malloc(L->n_archivos * 
					     sizeof(struct uvf_sample *));
  printf("ok2\n");
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
  printf ("Calculando sincos\n");
  calcularSinCos(L);
  printf ("L listo\n");
  return L;
}

/* Calcula la funcion L.
 * En pars los iesimos parametros son la posicion x, los i+1 la posicion y
 * y los i+2 las intensidades de los n centros de Voronoi.
 */
double L (funcL *fL, double *pars, double n){
  int i, j, k;
  double ret = 0, S = 0;
  double *im;

  fL->malla = newMallaVoronoi();
  //printf ("Ok L 1\n");
  im = malloc (fL->fg_image->npixels * sizeof (double));
  for (i = 0; i < n; i++) {
    insertarSitio (fL->malla, pars[3 * i], pars[3 * i + 1], pars[3 * i + 2]);
  }
  //printf ("Ok L 2\n");
  fL->imagen = toImage(fL->malla, fL->fg_image->size[0], fL->fg_image->size[1]);
  
  //printf ("Ok L 3\n");
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[0]; j++) {
      fL->fg_image->pixels[i + j * fL->fg_image->size[0]] = fL->imagen[i][j];
    }
  }
  //printf ("Ok L 4\n");
    for (k = 0; k < fL->n_archivos; k++) {
    mockcbi_sub(fL->cmb_image, fL->fg_image, 
		fL->header_obs[k], fL->samples_obs[k], fL->samples_mod[k], 
		fL->beam);
    for(i = 0; i < fL->header_obs[k]->nsamp; i++) {
      for(j = 0; j < fL->header_obs[k]->nchan; j++) {
	ret += ((SQR(fL->samples_mod[k][i].rdata[j * 3] 
		     - fL->samples_obs[k][i].rdata[j * 3])
		 + SQR(fL->samples_mod[k][i].rdata[j * 3 + 1] 
		       - fL->samples_obs[k][i].rdata[j * 3 + 1]))
		* fL->samples_mod[k][i].rdata[j *3 + 2]);
      }
    }
    }
  //printf ("Ok L 5\n");
  ret = ret / 2 - S;
  //printf ("Chi2 = %g\n", ret);
  for (i = 0; i < fL->fg_image->size[0]; i++) {
    free (fL->imagen[i]);
  }
  free (fL->imagen);
  free (im);
  eliminarMallaVoronoi(fL->malla);
  /*  for (k = 0; k < fL->n_archivos; k++) {
    free(fL->samples_mod[k]);
    }*/
  return ret;
}

void dL (funcL *L, double *pars, double *grad, double n){
  int i, j, k;
  double dVdI_R, dVdI_I;
  PoligonoVoronoi *pol;
  NodoLista *nodo;

  L->malla = newMallaVoronoi();
  for (i = 0; i < n; i++) {
    insertarSitio (L->malla, pars[3 * i], pars[3 * i + 1], pars[3 * i + 2]);
  }
  for (i = 0, nodo = L->malla->poligonos->primero; i < n; i++, nodo = nodo->sgte) {
    pol = (PoligonoVoronoi *) nodo->info;
    if ((pol->x != pars[3 * i]) || (pol->y != pars[3 * i + 1])) {
      fprintf (stderr, "ERROR: dL lista no coincide con arreglo. :P\n");
      exit (1);
    }
    //dVdI_R = 
  }
}

void calcularSinCos(funcL *L){
  int i, j, l, ix, iy;
  long k;    // k recorre xi
  double dx, dy, x, y;
  double lobs, mobs, obsra, obsdec, raimage, decimage, ra_offset = 0, cte;
  double arc;      // variables para el haz primario
  int x0, y0;            // pixel en el centro del haz primario.
  double kx1, kx2, kx3, kx4;

  raimage = L->fg_image->crval[0]*RPDEG;
  decimage = L->fg_image->crval[1]*RPDEG;

  L->sin = (double ****) malloc (L->n_archivos * sizeof (double ***));
  L->cos = (double ****) malloc (L->n_archivos * sizeof (double ***));
  for (i = 0; i < L->n_archivos; i++) {
    L->sin[i] = (double ***) malloc(L->fg_image->npixels * sizeof(double**));
    L->cos[i] = (double ***) malloc(L->fg_image->npixels * sizeof(double**));
    obsra = ra_offset + L->header_obs[i]->obsra*RPDEG;
    obsdec = L->header_obs[i]->obsdec*RPDEG;  
    /*if (debug)
      printf("aten: pointing center: %s, %s (%g arcmin from image center)\n",
	     ra_string(obsra), dec_string(obsdec),
	     slaDsep(raimage, decimage, obsra, obsdec)/RPARCM);*/
    //printf("obsra: %d, obsdec: %d\n", obsra, obsdec);
    direccos(obsra, obsdec, raimage, decimage, &lobs, &mobs);
    dx = L->fg_image->cdelt[0] * RPDEG;   // en radianes
    dy = L->fg_image->cdelt[1] * RPDEG;   // en radianes
    x0 = (L->fg_image->crpix[0]-1.0) + lobs/dx;
    y0 = (L->fg_image->crpix[1]-1.0) + mobs/dy;
    dx = L->fg_image->cdelt[0] * 60;   // en arcmin
    dy = L->fg_image->cdelt[1] * 60;   // en arcmin
    
    printf("crpix[0]: %d, crpix[1]: %d\n",
	   (int) L->fg_image->crpix[0] ,(int) L->fg_image->crpix[1]);
    printf("xo: %d, yo: %d\n", x0, y0);
    //  samples_mod = g_mockcbi(cmb_image, fg_image, header_obs, samples_obs);
    k = 1;
    for (iy = 0; iy < L->fg_image->size[1]; iy++) {
      y = (iy - y0) * dy;
      for (ix = 0; ix < L->fg_image->size[0]; ix++) {
	x = (ix - x0) * dx;
	arc = RPARCM * sqrt(x * x + y * y); // radio en radianes 
	L->sin[i][k] = (double **) malloc(L->header_obs[i]->nchan 
					  * sizeof(double*));
	L->cos[i][k] = (double **) malloc(L->header_obs[i]->nchan 
					  * sizeof(double*));
	
	for(j = 0; j < L->header_obs[i]->nchan; j++) {
	  L->sin[i][k][j] = (double *) malloc(L->header_obs[i]->nsamp 
					      * sizeof(double));
	  L->cos[i][k][j] = (double *) malloc(L->header_obs[i]->nsamp 
					       * sizeof(double));
	  for (l = 0; l < L->header_obs[i]->nsamp; l++) {
	    kx1 = (L->samples_mod[i][l].u * x 
		   + L->samples_mod[i][l].v * y) 
	      * L->header_obs[i]->iffreq[j];
	    kx2 = (L->samples_mod[i][l].u * (x + dx) 
		   + L->samples_mod[i][l].v * y) 
	      * L->header_obs[i]->iffreq[j];
	    kx3 = (L->samples_mod[i][l].u * x 
		   + L->samples_mod[i][l].v * (y + dy)) 
	      * L->header_obs[i]->iffreq[j];
	    kx4 = (L->samples_mod[i][l].u * (x + dx) 
		   + L->samples_mod[i][l].v * (y + dy)) 
	      * L->header_obs[i]->iffreq[j];

	    cte = -1.0 / (4.0 * PI * PI
			  * L->samples_mod[i][l].u * L->samples_mod[i][l].v
			  * L->header_obs[i]->iffreq[j] 
			  * L->header_obs[i]->iffreq[j]);
	    L->cos[i][k][j][l] = cte * (cos(2 * PI * kx1)
					- cos(2 * PI * kx2)
					- cos(2 * PI * kx3)
					+ cos(2 * PI * kx4));
	    L->sin[i][k][j][l] = cte * (sin(2 * PI * kx1)
					- sin(2 * PI * kx2)
					- sin(2 * PI * kx3)
					+ sin(2 * PI * kx4));
	    //seno = sin(2 * PI * kx);   
	    
	    //freq = header_obs[i]->iffreq[j] * 1e-9; 
	    //attenu[i][k][j] = primary_beam(arc, freq, &beam);
	  }
	}
	k++;
      }
    }
  }
}

void eliminarFuncL (funcL *L) {
  int i, j;
  Status status;

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
      free (L->atten[i][j]);
    }
    free (L->atten[i]);
  }
  free (L->atten);

  delete_map (L->cmb_image);
  delete_map (L->fg_image);
  
  free (L);
}
