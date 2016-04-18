#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include "frprmn.h"
#include "mcheck.h"
#include "reconstructor.h"

double gradienteNumerico (funcL *fL, double *pars, int i, double delta);
double *abrirArchivo (char *nombreArchivo, int *n);
void probarMalla(funcL *fL, double pars[]);

double gradienteNumerico (funcL *fL, double *pars, int i, double delta) {
  double fL1, fL2, par_old;
  //double *grad = (double *) malloc (3 * fL->n_pols * sizeof(double));

  par_old = pars[i];
  printf ("n gradienteNumerico = %d\n", fL->n_pols);
  fL1 = L (fL, pars, fL->n_pols);
  pars[i] += delta;
  fL2 = L (fL, pars, fL->n_pols);
  pars[i] = par_old;

  printf ("fl1 = %lf\n", fL1);
  printf ("fl2 = %lf\n", fL2);
  printf ("delta = %lf\n", delta);
  
  //return fL2 / delta - fL1 / delta;
  return (fL2 - fL1) / delta;
}

double *abrirArchivo (char *nombreArchivo, int *n) {
  double *pars;
  int i;
  FILE *archivo;

  archivo = fopen (nombreArchivo, "r");
  fscanf (archivo, "%d\n", n);
  printf ("n archivo = %d\n", *n);

  pars = (double *) malloc (3 * (*n) * sizeof(double));
  for (i = 0; i < *n; i++) {
    fscanf (archivo, "%lf\t%lf\t%lf\n", 
	    &pars[3 * i], &pars[3 * i + 1], &pars[3 * i + 2]);
    printf ("%lf\t%lf\t%lf\n", 
	    pars[3 * i], pars[3 * i + 1], pars[3 * i + 2]);
  }
  fclose (archivo);
  return pars;
}

int main (int argc, char **argv) {
  int n = 100, i, indice = 4;
  funcL *fL;
  FILE *out;
  double *pars, *grad, delta;
  Reconstructor *r;

  mtrace();
  
  if (argc < 2) {
    printf ("ERROR usage: comprobacion_gradiente archivo.in\n");
    exit (1);
  }
  
  r = leerArchivoEntrada (argv[1]);
  n = r->fL->n_pols;

  //pars = abrirArchivo (argv[argc - 1], &n);
  //pars = (double *) malloc (3 * n * sizeof(double));
  pars = (double *) calloc (3 * n, sizeof(double));
  srand (0);
  for (i = 0; i < n; i++) {
    pars[3 * i] = (double) random() / RAND_MAX;
    pars[3 * i + 1] = (double) random() / RAND_MAX;
    pars[3 * i + 2] = (double) random() / RAND_MAX;
    //pars[3 * i] = 3 * i / (3.0 * n);
    //pars[3 * i + 1] = (3 * i + 1)/ (3.0 * n);
    //pars[3 * i + 2] = (3 * i + 2)/ (3.0 * n);
  }

  //grad = (double *) malloc (3 * n * sizeof(double));
  grad = (double *) calloc (3 * n, sizeof(double));

  fL = r->fL;

  dL (fL, pars, grad, n, 1);

  do_write_fits(fL->fg_image, "!imagen.fits");
  
  out = fopen ("gradientex_exacto.dat", "w");

  for (i = 0; i < n; i++) {
    fprintf (out, "%d\t%.15g\n", i, grad[3 * i]);
  }
  
  fclose (out);

  out = fopen ("gradientey_exacto.dat", "w");

  for (i = 0; i < n; i++) {
    fprintf (out, "%d\t%.15g\n", i, grad[3 * i + 1]);
  }
  
  fclose (out);

  out = fopen ("gradienteI_exacto.dat", "w");

  for (i = 0; i < n; i++) {
    fprintf (out, "%d\t%.15g\n", i, grad[3 * i + 2]);
  }
  
  fclose (out);
  //exit(0);
  out = fopen ("gradientex.dat", "w");
  for (delta = 1; delta > 1e-6; delta /= 1.1) {
    double gradNum = gradienteNumerico(fL, pars, indice * 3, delta);
    fprintf (out, "%.16g\t%.16g\n", delta, gradNum);
    printf ("delta = %g,\tdL/dx_%d = %g\n", 
	     delta, indice, gradNum);
  }
  fclose (out);

  out = fopen ("gradientey.dat", "w");
  for (delta = 1; delta > 1e-6; delta /= 1.1) {
    double gradNum = gradienteNumerico(fL, pars, indice * 3 + 1, delta);
    fprintf (out, "%.16g\t%.16g\n", delta, gradNum);
    printf ("delta = %g,\tdL/dx_%d = %g\n", 
	     delta, indice, gradNum);
  }
  fclose (out);

  out = fopen ("gradienteI.dat", "w");
  for (delta = 1; delta > 1e-6; delta /= 1.1) {
    double gradNum = gradienteNumerico(fL, pars, indice * 3 + 2, delta);
    fprintf (out, "%.16g\t%.16g\n", delta, gradNum);
    printf ("delta = %g,\tdL/dx_%d = %g\n", 
	     delta, indice, gradNum);
  }
  fclose (out);
  
  for (indice = 0; indice < n; indice++) {
    char nombre[30];

    //eliminarFuncL (fL);
    //fL = newFuncL (argv[1], argv + 2, argc - 3, n);

    sprintf(nombre, "funcionx_%d.dat", indice);
    out = fopen (nombre, "w");
    for (delta = -0.01; delta < 0.01; delta += 1e-4) {
      double par_old = pars[indice * 3];
      double func;
      
      pars[indice * 3] -= delta;
      func = L (fL, pars, fL->n_pols);
      pars[indice * 3] = par_old;
      fprintf (out, "%.6g\t%.16g\n", delta, func);
      printf ("%g\t%.16g\n", delta, func);
    }
    fclose (out);  
    
    sprintf(nombre, "funciony_%d.dat", indice);
    out = fopen (nombre, "w");
    for (delta = -0.01; delta < 0.01; delta += 1e-4) {
      double par_old = pars[indice * 3 + 1];
      double func;
      
      pars[indice * 3 + 1] += delta;
      func = L (fL, pars, fL->n_pols);
      pars[indice * 3 + 1] = par_old;
      fprintf (out, "%.6g\t%.16g\n", delta, func);
      printf ("%g\t%.16g\n", delta, func);
    }
    fclose (out);  
    
    sprintf(nombre, "funcionI_%d.dat", indice);
    out = fopen (nombre, "w");
    for (delta = -0.01; delta < 0.01; delta += 1e-4) {
      double par_old = pars[indice * 3 + 2];
      double func;
      
      pars[indice * 3 + 2] += delta;
      func = L (fL, pars, fL->n_pols);
      pars[indice * 3 + 2] = par_old;
      fprintf (out, "%.6g\t%.16g\n", delta, func);
      printf ("%g\t%.16g\n", delta, func);
    }
    fclose (out);
  }
  
  eliminarFuncL (fL);
  free (pars);
  free (grad);
  
  return 0;
}


