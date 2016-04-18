#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include <gsl/gsl_multimin.h>

funcL *fL;
int iter_grad;

double my_f (const gsl_vector *v, void *params)
{
  double ret;
  int n = ((int *)params)[0], i;
  double *pars = malloc (n * sizeof(double));
  
  printf ("Calculando chi2\n");
  for (i = 0; i < n; i++) {
    /*
    if (gsl_vector_get(v, i) < 0) {
      gsl_vector_set(v, i, 0);
      }*/
    //pars[3 * i] = gsl_vector_get(v, 3 * i);
    //pars[3 * i + 1] = gsl_vector_get(v, 3 * i + 1);
    pars[i] = gsl_vector_get(v, i);
  }
  ret = L(fL, pars, n);
  printf ("chi2 = %.55g\n", ret);
  free (pars);
  return ret;
}

void my_df (const gsl_vector *v, void *params,
	    gsl_vector *df)
{
  int n = ((int *)params)[0], i;
  double *pars = malloc (n * sizeof(double)),
    *grad = malloc (n * sizeof(double));
  char nombre[20];
  FILE *archivo;
  
  sprintf (nombre, "grad_%d.dat", iter_grad);
  archivo = fopen (nombre, "w");

  printf ("n = %d\n", n);
  for (i = 0; i < n; i++) {
    /*
    if (gsl_vector_get(v, i) < 0) {
      gsl_vector_set(v, i, 0);
      }*/
    //pars[3 * i] = gsl_vector_get(v, 3 * i);
    //pars[3 * i + 1] = gsl_vector_get(v, 3 * i + 1);
    pars[i] = gsl_vector_get(v, i);
  }  
  printf ("Calculando gradiente\n");
  dL (fL, pars, grad, n);
  
  //printf ("grad: ");
  for (i = 0; i < n; i++) {
    //gsl_vector_set(v, 3 * i, grad[3 * i]);
    //gsl_vector_set(v, 3 * i + 1, grad[3 * i]);
    //printf (" %g", grad[i]);
    fprintf (archivo, "%d\t%.25g\n", i, grad[i]);
    gsl_vector_set(df, i, grad[i]);
  } 
  //printf ("\n");
  fclose (archivo);
}

void my_fdf (const gsl_vector *x, void *params,
	     double *f, gsl_vector *df)
{
  *f = my_f(x, params);
  my_df(x, params, df);
}

int minimizarI(char *fits_name, char **uvf_names, int n_uvfiles,
	       int n, double init_value, double step_size, double tol, double epsabs) {
  int par[1] = {n};
  size_t iter = 0, np = n, i;
  int status;
  
  char *fileout = (char*) malloc(30 * sizeof(char));

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;
  //double size;
  struct image *mask;
  FILE *archivoChi = fopen ("chi2.dat", "w");
  FILE *archivoCorr;
  char nombre[30];
  
  my_func.f = &my_f;
  my_func.df = &my_df;
  my_func.fdf = &my_fdf;
  my_func.n = np;
  my_func.params = &par;
  
  fL = newFuncL(fits_name, uvf_names, n_uvfiles, n);
  mask = do_read (fits_name);
  
  x = gsl_vector_alloc (np);
  for (i = 0; i < np; i++) {
    gsl_vector_set(x, i, init_value);
  }

  T = gsl_multimin_fdfminimizer_conjugate_pr;
  s = gsl_multimin_fdfminimizer_alloc (T, np);
  
  iter_grad = 0;
  gsl_multimin_fdfminimizer_set (s, &my_func, x, step_size, tol);
  fprintf (archivoChi, "%.25g\n", s->f);
  //printf ("%10.10e\n", s->f);
  fclose (archivoChi);
  do {
    NodoLista *nodo;
    PoligonoVoronoi *pol;

    sprintf (nombre, "corr_%d.dat", iter);
    archivoCorr = fopen (nombre, "w");
    for (i = 0, nodo = fL->malla->poligonos->primero; 
	 i < n; i++, nodo = nodo->sgte) {
      pol = (PoligonoVoronoi *) nodo->info;
      fprintf (archivoCorr, "%d\t%g\t%g\n", i, 
	       sqrt (pow(pol->x, 2) + pow(pol->y, 2)), pol->valor);
    }
    fclose (archivoCorr);

    iter_grad++;
    printf ("Iteracion: %d\n", iter);
    //if (iter % 1000 == 0) {
      sprintf(fileout, "!MEM_%d.fits", iter);
      do_write_fits(fL->fg_image, fileout);
      /*
      for (i = 0; i < mask->size[0]; i++) {
	for (j = 0; j < mask->size[0]; j++) {
	  mask->pixels[i + j * mask->size[0]] = fL->mask[i][j];
	  //printf ("mask[%d][%d] = %g\n", i, j, fL->mask[i][j]);
	}
      }
      sprintf(fileout, "!mask_%d.fits", iter);
      do_write_fits(mask, fileout);
    }*/
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);

    archivoChi = fopen ("chi2.dat", "a");
    fprintf (archivoChi, "%.25g\n", s->f);
    //printf ("%10.10e\n", s->f);
    fclose (archivoChi);
    
    if (status) {
      fprintf (stderr, "ERROR GSL: %s\n", gsl_strerror(status));
      break;
    }
    
    status = gsl_multimin_test_gradient (s->gradient, epsabs);
	
    if (status == GSL_SUCCESS) {
      printf ("Minimum found at:\n");
      }
  }
  while (status == GSL_CONTINUE && iter < 1e10);

  gsl_vector_free (x);
  gsl_multimin_fdfminimizer_free (s);
  
  return status;  
}

int main(int argc, char **argv) {
  int n = 500, status;

  // status = minimizarI (n, 0.0, 1e-2, 1e-4, 1e-3); // ROPHW
  status = minimizarI (argv[1], argv + 2, argc - 2, 
		       n, 1.0, 1e2, 1e-10, 1e-10);

    
  /*  for (i = 0; i < n; i++) {
      printf(" (%g, %g) = %g\n", gsl_vector_get(s->x, i),
      gsl_vector_get(s->x, 3*i + 1), gsl_vector_get(s->x, 3*i + 2));
      }*/
  do_write_fits(fL->fg_image, "!IC443_mem.fits");
  
  return status;
  
}
