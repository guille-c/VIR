#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arista.h"
#include "lista.h"
#include "mallaVoronoi.h"
#include "funcL.h"
#include <gsl/gsl_multimin.h>

funcL *f;


double my_f (const gsl_vector *v, void *params)
{
  double ret;
  int n = ((int *)params)[0], i;
  double *pars = malloc (3 * n * sizeof(double));
  
  for (i = 0; i < n; i++) {
    if (gsl_vector_get(v, 3 * i + 2) < 0) {
      gsl_vector_set(v, 3 * i + 2, 0);
    }
    pars[3 * i] = gsl_vector_get(v, 3 * i);
    pars[3 * i + 1] = gsl_vector_get(v, 3 * i + 1);
    pars[3 * i + 2] = gsl_vector_get(v, 3 * i + 2);
  }
  ret = L(f, pars, n);
  //ret = (double) random();
  free (pars);
  return ret;
}

int main(int argc, char **argv) {
  int n = 500, par[1] = {n};
  size_t iter = 0, np = n * 3, i, j;
  int status;
  
  char *fileout = (char*) malloc(30 * sizeof(char));
  const gsl_multimin_fminimizer_type *T = 
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *x, *ss;
  gsl_multimin_function my_func;
  double size;
  struct image *mask;
  
  //printf ("%d\n", argc);
  //printf ("NULL_PIX: %g\n", NULL_PIX);
  f = newFuncL(argv[1], argv + 2, argc - 2);
  mask = do_read (argv[1]);

  ss = gsl_vector_alloc (np);
  gsl_vector_set_all (ss, 1e-1);

  x = gsl_vector_alloc (np);
  for (i = 0; i < n; i++) {
    gsl_vector_set(x, 3 * i, (double) random() / RAND_MAX);
    gsl_vector_set(x, 3 * i + 1, (double) random() / RAND_MAX);
    gsl_vector_set(x, 3 * i + 2, 0.0);
  }
  my_func.f = &my_f;
  my_func.n = np;
  my_func.params = (void *) &par;

  s = gsl_multimin_fminimizer_alloc (T, np);
  gsl_multimin_fminimizer_set (s, &my_func, x, ss);
  
  do {
    printf ("Iteracion: %d\n", iter);
    if (iter % 1000 == 0) {
      sprintf(fileout, "!MEM_%d.fits", iter);
      do_write_fits(f->fg_image, fileout);
      for (i = 0; i < mask->size[0]; i++) {
	for (j = 0; j < mask->size[0]; j++) {
	  mask->pixels[i + j * mask->size[0]] = f->mask[i][j];
	  //printf ("mask[%d][%d] = %g\n", i, j, f->mask[i][j]);
	}
      }
      sprintf(fileout, "!mask_%d.fits", iter);
      do_write_fits(mask, fileout);
    }
    iter++;
    status = gsl_multimin_fminimizer_iterate (s);
    
    if (status)
      break;
    
    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, 1e-3);
	
    if (status == GSL_SUCCESS) {
      printf ("Minimum found at:\n");
      }
  }
  while (status == GSL_CONTINUE && iter < 1e10);
  
  for (i = 0; i < n; i++) {
    printf(" (%g, %g) = %g\n", gsl_vector_get(s->x, 3*i),
	   gsl_vector_get(s->x, 3*i + 1), gsl_vector_get(s->x, 3*i + 2));
  }
  gsl_vector_free (x);
  gsl_vector_free (ss);
  gsl_multimin_fminimizer_free (s);
  
  do_write_fits(f->fg_image, "!IC443_mem.fits");
  
  return status;

}
