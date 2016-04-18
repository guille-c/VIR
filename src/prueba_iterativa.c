#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include "frprmn.h"
#include <gsl/gsl_multimin.h>

funcL *fL;
int iter_grad;

//----------I NR ROUTINES-------------------------------------------

float chi2NumRecI (float pars[]) {
  double ret; 
  double *parsD = (double *) malloc (fL->n_pols * sizeof (double));
  int i;

  for (i = 1; i <= fL->n_pols; i++) {
    parsD [i - 1] = (double) pars[i];
  }
  
  printf ("Calculando chi2\n");
  ret  =  L(fL, parsD, fL->n_pols);
  printf ("chi2 = %.55g\n", ret);
  
  return ret;
}

void dchi2NumRecI (float pars[], float grad[]) {
  double *parsD = (double *) malloc (fL->n_pols * sizeof (double));
  double *gradD = (double *) malloc (fL->n_pols * sizeof (double));
  int i;
  char fileout[30];

  iter_grad ++;
  for (i = 1; i <= fL->n_pols; i++) {
    parsD [i - 1] = (double) pars[i];
  }
  printf ("Calculando dchi2\n");
  dL(fL, parsD, gradD, fL->n_pols);

  for (i = 1; i <= fL->n_pols; i++) {
    grad [i] = (float) gradD [i - 1];
  }

  sprintf(fileout, "!MEM_%d.fits", iter_grad);
  do_write_fits(fL->fg_image, fileout);

  printf ("dchi2 calculado\n");
}

double minimizarI_NR (char *fits_name, char **uvf_names, int n_uvfiles,
		      int n, double init_value,
		      float ftol) {
  int iter, i;
  float fret;
  struct image *mask;
  long t;
  float *p = (float *) malloc (n * sizeof (float));

  iter_grad = 0;
  for (i = 0; i < n; i++) {
    p[i] = init_value;
    //printf ("p[%d] = %f\n", i, p[i]);
  }

  fL = newFuncL(fits_name, uvf_names, n_uvfiles, n);
  mask = do_read (fits_name);

  t = time(0);
  frprmn(p - 1, n, ftol, &iter, &fret, chi2NumRecI, dchi2NumRecI);
  t = time(0) - t;
  printf("Terminado P-R en %d segundos y %d iteraciones\n", (int) t, iter);

  do_write_fits(fL->fg_image, "!MEM_NR.fits");
  return fret;
}

//----------(x, y) NR ROUTINES---------------------------------------------

float chi2NumRec (float pars[]) {
  double ret; 
  double *parsD = (double *) malloc (3 * fL->n_pols * sizeof (double));
  int i;

  for (i = 0; i < fL->n_pols; i++) {
    parsD [3 * i]     = (double) pars[3 * i + 1];
    parsD [3 * i + 1] = (double) pars[3 * i + 1 + 1];
    parsD [3 * i + 2] = (double) pars[3 * i + 2 + 1];
  }
  
  printf ("Calculando chi2\n");
  ret  =  L(fL, parsD, fL->n_pols);
  printf ("chi2 = %.55g\n", ret);
  
  return ret;
}

void dchi2NumRec (float pars[], float grad[]) {
  double *parsD = (double *) malloc (3 * fL->n_pols * sizeof (double));
  double *gradD = (double *) malloc (3 * fL->n_pols * sizeof (double));
  int i, j;
  char fileout[30];
  struct image *imagenPix;

  iter_grad ++;
  for (i = 0; i < fL->n_pols; i++) {
    parsD [3 * i]     = (double) pars[3 * i + 1];
    parsD [3 * i + 1] = (double) pars[3 * i + 1 + 1];
    parsD [3 * i + 2] = (double) pars[3 * i + 2 + 1];
    /*
    printf ("parsD[%d] = %g\n", 3 * i, parsD[3 * i]);
    printf ("parsD[%d] = %g\n", 3 * i + 1, parsD[3 * i + 1]);
    printf ("parsD[%d] = %g\n\n", 3 * i + 2, parsD[3 * i + 2]);
    */
  }
  printf ("Calculando dchi2\n");
  dL(fL, parsD, gradD, fL->n_pols);

  for (i = 0; i < fL->n_pols; i++) {
    grad [3 * i + 1] = (float) gradD [3 * i];
    grad [3 * i + 1 + 1] = (float) gradD [3 * i + 1];
    grad [3 * i + 2 + 1] = (float) gradD [3 * i + 2];
    //printf ("grad[%d] = %g, grad[%d] = %g, grad[%d] = %g\n",
    //    3 * i + 1, grad [3 * i + 1], 
    //    3 * i + 1 + 1, grad [3 * i + 1 + 1], 
    //    3 * i + 1 + 2, grad [3 * i + 1 + 2]);
  }

  sprintf(fileout, "!MEM_%d.fits", iter_grad);
  do_write_fits(fL->fg_image, fileout);

  for (i = 0; i < fL->fg_image->size[0]; i++) {
    for (j = 0; j < fL->fg_image->size[1]; j++) {
      imagenPix->pixels[i + j * fL->fg_image->size[0]] = fL->PixIntegrados_x[i][j];
    }
  }

  imagenPix = do_read(fL->nombreFits);
  sprintf (fileout, "!PixIntegrados_x%d.fits", iter_grad);
  do_write_fits(imagenPix, fileout);
  delete_map (imagenPix);

  printf ("dchi2 calculado\n");
}

double minimizar_NR (char *fits_name, char **uvf_names, int n_uvfiles,
		     int n, double init_value,
		     float ftol) {
  int iter, i;
  float fret;
  struct image *mask;
  long t;
  float *p = (float *) malloc (3 * n * sizeof (float));

  iter_grad = 0;
  for (i = 0; i < n; i++) {
    p[3 * i]     = (float) random() / RAND_MAX;
    p[3 * i + 1] = (float) random() / RAND_MAX;
    p[3 * i + 2] = init_value;
  }

  fL = newFuncL(fits_name, uvf_names, n_uvfiles, n);
  mask = do_read (fits_name);

  t = time(0);
  frprmn(p - 1, 3 * n, ftol, &iter, &fret, chi2NumRec, dchi2NumRec);
  t = time(0) - t;
  printf("Terminado P-R en %d segundos y %d iteraciones\n", (int) t, iter);

  do_write_fits(fL->fg_image, "!MEM_NR.fits");
  return fret;
}

int main(int argc, char **argv) {
  int n = 200;
  float ftol = 1e-10;
 
  minimizar_NR (argv[1], argv + 2, argc - 2, n, 0.0, ftol);
  return 1;
}
