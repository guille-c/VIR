/*********************
 * Comprobacion dVdx *
 *********************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include "frprmn.h"
#include "reconstructor.h"
#include "mcheck.h"

double gradienteNumerico (funcL *fL, double *pars, int i, double delta);
double *abrirArchivo (char *nombreArchivo, int *n);
void probarMalla(funcL *fL, double pars[]);

int main (int argc, char **argv) {
  int n = 10, i, indice = 4;
  funcL *fL;
  FILE *outxR, *outxI, *outyR, *outyI;
  double *pars, *grad, delta;
  double *dVdx_R, *dVdx_I, *dVdy_R, *dVdy_I, 
    datamin = 1e300, datamax = -1e300;
  int arch = 0, samp = 1, iff = 1;
  NodoLista *nodo;
  PoligonoVoronoi *pol;
  int entropia = 0;
  Reconstructor *r;
  //int aprox = 1;
  //mtrace();
  
  if (argc < 2) {
    printf ("ERROR usage: comprobacion_gradiente archivo.in\n");
    exit (1);
  }
  
  r = leerArchivoEntrada (argv[1]);
  n = r->fL->n_pols;
  //pars = abrirArchivo (argv[argc - 1], &n);
  //pars = (double *) malloc (3 * n * sizeof(double));
  dVdx_R = (double *) malloc (3 * n * sizeof(double));
  dVdx_I = (double *) malloc (3 * n * sizeof(double));
  dVdy_R = (double *) malloc (3 * n * sizeof(double));
  dVdy_I = (double *) malloc (3 * n * sizeof(double));
  pars = (double *) calloc (3 * n, sizeof(double));
  srand (0);
  for (i = 0; i < n; i++) {
    pars[3 * i] = (double) random() / RAND_MAX;
    //pars[3 * i] = 1.0/(2*n) + i * 1.0/n;
    //pars[3 * i] = 0.25 + i * 0.5;
    pars[3 * i + 1] = (double) random() / RAND_MAX;
    //pars[3 * i + 1] = 0.5 + i *0.1;
    //pars[3 * i + 1] = 0.5;
    pars[3 * i + 2] = (double) random() / RAND_MAX;
    //pars[3 * i + 2] = 1.0;
    //pars[3 * i] = 3 * i / (3.0 * n);
    //pars[3 * i + 1] = (3 * i + 1)/ (3.0 * n);
    //pars[3 * i + 2] = (3 * i + 2)/ (3.0 * n);
    printf ("pars%d = (%g, %g, %g)\n", 
	    i, pars[3 * i], pars[3 * i + 1], pars[3 * i + 2]);
  }

  //grad = (double *) malloc (3 * n * sizeof(double));
  grad = (double *) calloc (3 * n, sizeof(double));

  //fL = newFuncL (argv[1], argv + 2, argc - 2, n, entropia);
  fL = r->fL;
  L (fL, pars, fL->n_pols);  
  printf ("1) u, v = (%g, %g)\n", 
	  fL->samples_mod[arch][samp].u, fL->samples_mod[arch][samp].v);
  guardarFits (fL->imagen, fL->mask, "dV", fL->nombreFits);
  //exit (0);

  outxR = fopen ("dVdxR_exacto.dat", "w");
  outxI = fopen ("dVdxI_exacto.dat", "w");
  outyR = fopen ("dVdyR_exacto.dat", "w");
  outyI = fopen ("dVdyI_exacto.dat", "w");

  for (i = 0, nodo = fL->malla->poligonos->primero; 
       i < n; i++, nodo = nodo->sgte) {
    pol = (PoligonoVoronoi *) nodo->info;
    if (n - 1 - i != pol->id - 3) {
      fprintf (stderr, "ERROR: pol->id = %d != %d\n", 
	       pol->id - 3, n - 1 - i);
      exit (1);
    }
    dVdx (fL, arch, iff, samp, pol, 
	  &(dVdx_R[pol->id - 3]), &(dVdx_I[pol->id - 3]), 
	  &(dVdy_R[pol->id - 3]), &(dVdy_I[pol->id - 3]), 0);
  }
  for (i = 0; i < n; i++) {
    fprintf (outxR, "%d\t%g\n", i, dVdx_R[i]);
    fprintf (outxI, "%d\t%g\n", i, dVdx_I[i]);
    fprintf (outyR, "%d\t%g\n", i, dVdy_R[i]);
    fprintf (outyI, "%d\t%g\n", i, dVdy_I[i]);
  }
  fclose (outxR);
  fclose (outxI);
  fclose (outyR);
  fclose (outyI);

  outxR = fopen ("dVdxR_aprox.dat", "w");
  outxI = fopen ("dVdxI_aprox.dat", "w");
  outyR = fopen ("dVdyR_aprox.dat", "w");
  outyI = fopen ("dVdyI_aprox.dat", "w");

  for (i = 0, nodo = fL->malla->poligonos->primero; 
       i < n; i++, nodo = nodo->sgte) {
    pol = (PoligonoVoronoi *) nodo->info;
    if (n - 1 - i != pol->id - 3) {
      fprintf (stderr, "ERROR: pol->id = %d != %d\n", 
	       pol->id - 3, n - 1 - i);
      exit (1);
    }
    dVdx (fL, arch, iff, samp, pol, 
	  &(dVdx_R[pol->id - 3]), &(dVdx_I[pol->id - 3]), 
	  &(dVdy_R[pol->id - 3]), &(dVdy_I[pol->id - 3]), 1);
  }
  for (i = 0; i < n; i++) {
    fprintf (outxR, "%d\t%g\n", i, dVdx_R[i]);
    fprintf (outxI, "%d\t%g\n", i, dVdx_I[i]);
    fprintf (outyR, "%d\t%g\n", i, dVdy_R[i]);
    fprintf (outyI, "%d\t%g\n", i, dVdy_I[i]);
  }
  fclose (outxR);
  fclose (outxI);
  fclose (outyR);
  fclose (outyI);
  //exit(0);
  
  for (indice = 0; indice < n; indice++) {
    char nombre[30];
    int cont = 0;

    //eliminarFuncL (fL);
    //fL = newFuncL (argv[1], argv + 2, argc - 3, n);

    sprintf(nombre, "VxR_%d.dat", indice);
    outxR = fopen (nombre, "w");
    sprintf(nombre, "VxI_%d.dat", indice);
    outxI = fopen (nombre, "w");

    sprintf(nombre, "VyR_%d.dat", indice);
    outyR = fopen (nombre, "w");
    sprintf(nombre, "VyI_%d.dat", indice);
    outyI = fopen (nombre, "w");

    for (delta = -0.1; delta < 0.1; delta += 1e-3) {
      double par_old = pars[indice * 3];
      double func;

      cont++;

      pars[indice * 3] -= delta;
      func = L (fL, pars, fL->n_pols);
      
      if (cont % 10 == 0) {
	if (cont < 10) {
	  sprintf(nombre, "mallax_%d_00%d", indice, cont);
	}
	else if (cont < 100) {
	  sprintf(nombre, "mallax_%d_0%d", indice, cont);
	}
	else {
	  sprintf(nombre, "mallax_%d_%d", indice, cont);
	}
	guardarFits (fL->imagen, fL->mask, nombre, fL->nombreFits);
      }
      pars[indice * 3] = par_old;
      fprintf (outxR, "%.6g\t%.16g\n", delta, 
	       fL->samples_mod[arch][samp].rdata[iff * 3]);
      fprintf (outxI, "%.6g\t%.16g\n", delta, 
	       fL->samples_mod[arch][samp].rdata[iff * 3 + 1]);
      //printf ("%g\t%.16g\n", delta, func);

      par_old = pars[indice * 3 + 1];
      pars[indice * 3 + 1] += delta;
      func = L (fL, pars, fL->n_pols);
      
      if (cont % 10 == 0) {
	if (cont < 10) {
	  sprintf(nombre, "mallay_%d_00%d", indice, cont);
	}
	else if (cont < 100) {
	  sprintf(nombre, "mallay_%d_0%d", indice, cont);
	}
	else {
	  sprintf(nombre, "mallay_%d_%d", indice, cont);
	}
	guardarFits (fL->imagen, fL->mask, nombre, fL->nombreFits);
      }

      pars[indice * 3 + 1] = par_old;
      fprintf (outyR, "%.6g\t%.16g\n", delta, 
	       fL->samples_mod[arch][samp].rdata[iff * 3]);
      fprintf (outyI, "%.6g\t%.16g\n", delta, 
	       fL->samples_mod[arch][samp].rdata[iff * 3 + 1]);
      //printf ("%g\t%.16g\n", delta, func);
    }
    fclose (outxR);
    fclose (outxI);
    fclose (outyR);
    fclose (outyI);
  }
  
  eliminarFuncL (fL);
  free (pars);
  free (grad);
  
  return 0;
}


