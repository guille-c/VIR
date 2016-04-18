/*********************
 * Comprobacion dVdI *
 *********************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include "frprmn.h"
#include "funcionBayesVoronoiCBI.h"
#include "mcheck.h"

double gradienteNumerico (funcL *fL, double *pars, int i, double delta);
double *abrirArchivo (char *nombreArchivo, int *n);
void probarMalla(funcL *fL, double pars[]);
void compararVisibilidades(funcL *fL, double *pars);
void imprimirVisibilidades (funcL *fL, char *nombre_archivo);

int main (int argc, char **argv) {
  int n = 10, i, j, indice = 4, idPol;
  funcL *fL;
  FILE *outR, *outI;
  double *pars, *grad, delta;
  double *dVdI_R, *dVdI_I,
    datamin = 1e300, datamax = -1e300;
  int arch = 0, samp = 1, iff = 1;
  NodoLista *nodo;
  PoligonoVoronoi *pol;
  int entropia = 0;
  double u, v, x, y, kx;
  struct image *mask10;
  FuncionBayesVoronoiCBI *func = new FuncionBayesVoronoiCBI (argv[1]);
  //int aprox = 1;
  //mtrace();
  
  if (argc < 2) {
    printf ("ERROR usage: comprobacion_gradiente archivo.in\n");
    exit (1);
  }
  
  fL = func->getFL();

  mask10 = do_read(fL->nombreFits);
  for (i = 0; i < mask10->npixels; i++) {
    mask10->pixels[i] = 0;
  }

  n = fL->n_pols;
  //pars = abrirArchivo (argv[argc - 1], &n);
  //pars = (double *) malloc (3 * n * sizeof(double));
  dVdI_R = (double *) malloc (n * sizeof(double));
  dVdI_I = (double *) malloc (n * sizeof(double));
  pars = (double *) calloc (3 * n, sizeof(double));
  srand (0);
  for (i = 0; i < n; i++) {
    dVdI_R[i] = 0;
    dVdI_I[i] = 0;
    pars[3 * i] = (double) random() / RAND_MAX;
    //pars[3 * i] = 1.0/(2*n) + i * 1.0/n;
    //pars[3 * i] = 0.25 + i * 0.5;
    pars[3 * i + 1] = (double) random() / RAND_MAX;
    //pars[3 * i + 1] = 0.5 + i *0.1;
    //pars[3 * i + 1] = 0.5;
    pars[3 * i + 2] = (double) random() / RAND_MAX * 10.5;
    //pars[3 * i + 2] = 1.0;
    //pars[3 * i] = 3 * i / (3.0 * n);
    //pars[3 * i + 1] = (3 * i + 1)/ (3.0 * n);
    //pars[3 * i + 2] = (3 * i + 2)/ (3.0 * n);
    printf ("pars%d = (%g, %g, %g)\n", 
	    i, pars[3 * i], pars[3 * i + 1], pars[3 * i + 2]);
  }

  compararVisibilidades(fL, pars);
  exit(0);

  //grad = (double *) malloc (3 * n * sizeof(double));
  grad = (double *) calloc (3 * n, sizeof(double));

  //fL = newFuncL (argv[1], argv + 2, argc - 2, n, entropia);
  L (fL, pars, fL->n_pols, MOCKCBI);
  printf ("1) u, v = (%g, %g)\n", 
	  fL->samples_mod[arch][samp].u, fL->samples_mod[arch][samp].v);
  guardarFits (fL->imagen, fL->mask, "dV", fL->nombreFits);
  //exit (0);

  dL (fL, pars, grad, fL->n_pols, EXACTA, MOCKCBI);

  outR = fopen ("dVdIR_exacto.dat", "w");
  outI = fopen ("dVdII_exacto.dat", "w");

  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      if (fL->mask[i][j] < 0) {
	fprintf (stderr, "ERROR en dLdI: fL->mask[%d][%d] = %d\n", 
		 i, j, fL->mask[i][j]);
	exit (1);
      }
      y = (j - fL->y0[arch]) * fL->dy[arch];
      x = (i - fL->x0[arch]) * fL->dx[arch];
      //printf ("x = %g, y = %g, %g", x, y, sqrt(1 - x*x - y*y)); exit(0);
      u = fL->samples_mod[arch][samp].u * fL->header_obs[arch]->iffreq[iff];
      v = fL->samples_mod[arch][samp].v * fL->header_obs[arch]->iffreq[iff];
      kx = (u * x + v * y);

      idPol = fL->mask[i][j] - 3;
      if (idPol == 10) {
	//printf ("sumando (%d, %d)\n", i + 1, j + 1);
	mask10->pixels[i + j * fL->fg_image->size[0]] ++;
      }
      //grad[idPol] += fL->fourierI[i][j] * fL->difmapNoise * fL->fg_scale;
      //printf("dVdI_R[%d] = %g\n", idPol, dVdI_R[idPol]);
      //dVdI_R[idPol] += fL->sin[arch][samp][iff] 
      dVdI_R[idPol] += fabs(fL->dx[arch] * fL->dy[arch])
	* fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff] * cos (2 * PI * kx)
	* fL->difmapNoise * fL->fg_scale;

      //dVdI_I[idPol] += fL->sin[arch][samp][iff] 
      dVdI_I[idPol] += fabs(fL->dx[arch] * fL->dy[arch])
	* fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff] * sin (2 * PI * kx)
	* fL->difmapNoise * fL->fg_scale; 

//      printf ("four = %g, atten = %g, cos = %g, dVdIR[%d] = %g\n", fL->sin[arch][samp][iff],
//	      fL->atten[arch][i + j * fL->fg_image->size[0] + 1][iff], cos (2 * PI * kx), 
//	      idPol, dVdI_R[idPol]);
      //grad[3 * idPol + 2] += fL->fourierI[i][j] * fL->difmapNoise * fL->fg_scale * areas[idPol] / areasPix[idPol]; Esto NO
      //grad[3 * idPol + 2] += fL->fourierI[i][j] * fL->difmapNoise * fL->fg_scale / areas[idPol] ;
    }
  }

  do_write_fits (mask10, "!mask10.fits");

  for (i = 0; i < n; i++) {
    fprintf (outR, "%d\t%g\n", i, dVdI_R[i]);
    fprintf (outI, "%d\t%g\n", i, dVdI_I[i]);
  }
  fclose (outR);
  fclose (outI);
  
//   for (indice = 0; indice < n; indice++) {
//     char nombre[30];
//     int cont = 0;
// 
//     //eliminarFuncL (fL);
//     //fL = newFuncL (argv[1], argv + 2, argc - 3, n);
// 
//     sprintf(nombre, "VR_%d.dat", indice);
//     outR = fopen (nombre, "w");
//     sprintf(nombre, "VI_%d.dat", indice);
//     outI = fopen (nombre, "w");
// 
//     for (delta = -1; delta < 1; delta += 1e-2) {
//       double par_old = pars[indice * 3 + 2];
//       double func;
// 
//       pars[indice * 3 + 2] += delta;
//       func = L (fL, pars, fL->n_pols);
//       
// //      if (cont % 10 == 0) {
// //	if (cont < 10) {
// //	  sprintf(nombre, "mallaI_%d_00%d", indice, cont);
// //	}
// //	else if (cont < 100) {
// //	  sprintf(nombre, "mallaI_%d_0%d", indice, cont);
// //	}
// //	else {
// //	  sprintf(nombre, "mallaI_%d_%d", indice, cont);
// //	}
// //	guardarFits (fL->imagen, fL->mask, nombre, fL->nombreFits);
// //      }
// 
//       pars[indice * 3 + 2] = par_old;
//       fprintf (outR, "%.6g\t%.16g\n", delta, 
// 	       fL->samples_mod[arch][samp].rdata[iff * 3]);
//       fprintf (outI, "%.6g\t%.16g\n", delta, 
// 	       fL->samples_mod[arch][samp].rdata[iff * 3 + 1]);
//       //printf ("%g\t%.16g\n", delta, func);
//     }
//     fclose (outR);
//     fclose (outI);
//   }

  for (indice = 0; indice < n; indice++) {
    char nombre[30];
    int cont = 0;

    printf ("Indice = %d\n", indice);
    //eliminarFuncL (fL);
    //fL = newFuncL (argv[1], argv + 2, argc - 3, n);

    sprintf(nombre, "VR_2_%d.dat", indice);
    outR = fopen (nombre, "w");
    sprintf(nombre, "VI_2_%d.dat", indice);
    outI = fopen (nombre, "w");

    for (delta = -1; delta < 1; delta += 1e-2) {
      double par_old = pars[indice * 3 + 2];
      double func, visR, visI;

      pars[indice * 3 + 2] += delta;
      func = L (fL, pars, fL->n_pols, MOCKCBI);
      
      pars[indice * 3 + 2] = par_old;
      visMod (fL, arch, iff, samp, &visR, &visI);
      fprintf (outR, "%.6g\t%.16g\n", delta, visR);
      fprintf (outI, "%.6g\t%.16g\n", delta, visI);
      //printf ("%.6g\t%.16g\n", delta, visR);
      //printf ("%g\t%.16g\n", delta, func);
    }
    fclose (outR);
    fclose (outI);
  }

  
  eliminarFuncL (fL);
  free (pars);
  free (grad);
  
  return 0;
}

void compararVisibilidades(funcL *fL, double *pars) {
  double *visR =  new double [fL->n_archivos * fL->header_obs[0]->nif * 
			     fL->header_obs[0]->nsamp * sizeof(double)];
  double *visI =  new double [fL->n_archivos * fL->header_obs[0]->nif * 
			     fL->header_obs[0]->nsamp * sizeof(double)];
  int i = 0;
  FILE *archivo = fopen ("comparacionVis.dat", "w");
  FILE *archivoR = fopen ("comparacionVisR.dat", "w");
  FILE *archivoI = fopen ("comparacionVisI.dat", "w");
  
  L (fL, pars, fL->n_pols, MOCKCBI);
  imprimirVisibilidades (fL, "vis_mockcbi.dat");
  for (int arch = 0; arch < fL->n_archivos; arch++) {
    for (int iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
      for (int samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	if (fL->samples_mod[arch][samp].rdata[3 * iff + 1] != 0) {
	  visR[i] = fL->samples_mod[arch][samp].rdata[3 * iff];
	  visI[i] = fL->samples_mod[arch][samp].rdata[3 * iff + 1];
	}
	i++;
      }
    }
  }
  cout << "i = " << i << ", " << fL->n_archivos * fL->header_obs[0]->nif * 
    fL->header_obs[0]->nsamp << "\n";

  L (fL, pars, fL->n_pols, MOCK_EXACTO);
  imprimirVisibilidades (fL, "vis_exacto.dat");
  i = 0;
  for (int arch = 0; arch < fL->n_archivos; arch++) {
    for (int iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
      for (int samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	if (fL->samples_mod[arch][samp].rdata[3 * iff + 1] != 0) {
	  fprintf (archivo, "%g\t%g\n", sqrt (visR[i] * visR[i] + visI[i] * visI[i]),
		   sqrt (fL->samples_mod[arch][samp].rdata[3 * iff] * 
			 fL->samples_mod[arch][samp].rdata[3 * iff] + 
			 fL->samples_mod[arch][samp].rdata[3 * iff + 1] *
			 fL->samples_mod[arch][samp].rdata[3 * iff + 1]));
	  fprintf (archivoR, "%g\t%g\n", visR[i], fL->samples_mod[arch][samp].rdata[3 * iff]);
	  fprintf (archivoI, "%g\t%g\n", visI[i], fL->samples_mod[arch][samp].rdata[3 * iff + 1]);
	  
	}
	i++;
      }
    }
  }
  fclose (archivo);
}

void imprimirVisibilidades (funcL *fL, char *nombre_archivo){
  double vis, radio;
  FILE *archivo = fopen (nombre_archivo, "w");

  for (int arch = 0; arch < fL->n_archivos; arch++) {
    for (int iff = 0; iff < fL->header_obs[arch]->nif; iff++) {
      for (int samp = 0; samp < fL->header_obs[arch]->nsamp; samp++) {
	if (fL->samples_mod[arch][samp].rdata[3 * iff + 1] != 0) {
	  radio = sqrt(fL->samples_mod[arch][samp].u * fL->samples_mod[arch][samp].u +
		       fL->samples_mod[arch][samp].v * fL->samples_mod[arch][samp].v)
	    * fL->header_obs[arch]->iffreq[iff];
	  vis = sqrt (fL->samples_mod[arch][samp].rdata[3 * iff] * 
		      fL->samples_mod[arch][samp].rdata[3 * iff] + 
		      fL->samples_mod[arch][samp].rdata[3 * iff + 1] *
		      fL->samples_mod[arch][samp].rdata[3 * iff + 1]);
	  fprintf (archivo, "%g\t%g\n", radio, vis);
	}
	else {
	  printf ("w = 0\n");
	}
      }
    }
  }
  
  fclose (archivo);
}
