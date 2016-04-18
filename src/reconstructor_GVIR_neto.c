#include "reconstructor_GVIR_neto.h"

float func_max;
int iter_pikaia;
Reconstructor *r_out;
// ---------Funciones para pikaia-------------------------------
extern void pikaia_(float *(*ff)(int *n, float x[*n]), int *n, float ctrl[12], 
	 	        float x[*n], float *f, int *status);

extern void rninit_(int *seed);

float ret;
float *ff_(int *n, float x[*n]){
  int i;
  char fileout [30]; 
  double *parsD = (double *) malloc ((*n) * sizeof (double));
  
  imprimirLog ("reconstructor.log", 
	       "        Calculando funcion para la %d iteracion con parametros:\n",
	       r_out->iter);
  for (i = 0; i < r_out->fL->n_pols; i++) {
    parsD[3 * i]     = (double) x[3 * i];
    parsD[3 * i + 1] = (double) x[3 * i + 1];
    parsD[3 * i + 2] = (double) x[3 * i + 2] * r_out->Nmax;
    imprimirLog ("reconstructor.log", "        %g\t%g\t%g\n", 
		 parsD[3 * i], parsD[3 * i + 1], parsD[3 * i + 2]);
  }
  
  ret  =  -L (r_out->fL, parsD, r_out->fL->n_pols);

  imprimirLog ("funcion.log", "%d\t%g\n", iter_pikaia, ret);
  if (ret > func_max){
    sprintf(fileout, "pikaia_%d", iter_pikaia);
    guardarFits (r_out->fL->imagen, r_out->fL->mask, 
		 fileout, r_out->fL->nombreFits);
    imprimirLog ("pikaia.log", "%d\t%g\n", iter_pikaia, ret);
    func_max = ret;
    iter_pikaia++;
  }
  
  free (parsD);
  return &ret;
}
//---------------------------------------------------------
  
Reconstructor *newReconstructor (char * nombreImagen, 
				 char **nombresVis, int nVis, int n,
				 double init_value, int init_gauss,
				 int entropia, double cuantaSize, 
				 double cutoff) {
  int i, j, k, nx, ny, imax = 0, jmax = 0;
  FILE* archivoLog = fopen ("reconstructor.log", "w");
  fclose (archivoLog);
  archivoLog = fopen ("pikaia.log", "w");
  fclose (archivoLog);
  archivoLog = fopen ("funcion.log", "w");
  fclose (archivoLog);

  Reconstructor *r = (Reconstructor *) malloc (sizeof (Reconstructor));

  func_max = -1e15;
  iter_pikaia = 0;
  r->entropia = entropia;
  r->nVis = nVis;
  printf ("nVis = %d\n", nVis);
  printf ("r->nVis = %d\n\n", r->nVis);
  r->iter = 0;
  r->Acutoff = cutoff;
  r->Nmax = 1;
  r->individuos = 100;
  r->generaciones = 500;
  r->fL = newFuncL (nombreImagen, nombresVis, nVis, n, r->entropia, cuantaSize);
  r->nombreImagen = (char *) malloc (256 * sizeof (char));
  strcpy(r->nombreImagen, nombreImagen);
  r->nombresVis = (char **) malloc (nVis * sizeof (char *));
  imprimirLog ("reconstructor.log", "Reconstruyendo ");
  for (i = 0; i < nVis; i++){
    r->nombresVis[i] = (char *) malloc (256 * sizeof (char));
    strcpy(r->nombresVis[i], nombresVis[i]);
    imprimirLog ("reconstructor.log", "%s ", r->nombresVis[i]);
  }
  // Inicializamos la Imagen

  r->p = (float *) malloc (3 * r->fL->n_pols * sizeof (float));
  r->iter = 0;

  nx = r->fL->fg_image->size[0];
  ny = r->fL->fg_image->size[1];

  imprimirLog ("reconstructor.log", "\nutilizando como imagen de referencia%s\n ",
	       r->nombreImagen);
  imprimirLog ("reconstructor.log", "nx = %d, ny = %d\n ", nx, ny);
  imprimirLog ("reconstructor.log", "entropia: %d\n", r->entropia);
  if (init_gauss) {
    double mediax, sigmax, mediay, sigmay, atten_max = -1e300;
    Gaussiana *gx, *gy;
    imprimirLog ("reconstructor.log", "inicializando con distribucion gaussiana\n");

    // Buscamos el maximo en atten.
    k = 1;
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
	if (atten_max < r->fL->atten[0][k][0]){
	  atten_max = r->fL->atten[0][k][0];
	  imax = i;
	  jmax = j;
	}
	k++;
      }
    }

    mediax = (double) imax / (nx - 1);
    mediay = (double) jmax / (ny - 1);
    
    for (i = imax; i < nx && 
	   r->fL->atten[0][jmax * nx + i + 1][0] > atten_max / 2; i++);
    sigmax = sqrt((double) (i - imax) / (nx - 1)); // Problema: se
						   // salen de la
						   // imagen
    if (sigmax > 0.2) {
      sigmax = 0.2;
    }
    
    for (j = jmax; j < ny && 
	   r->fL->atten[0][j * nx + imax + 1][0] > atten_max / 2; j++);
    sigmay = sqrt((double) (j - jmax) / (ny - 1));
    if (sigmay > 0.2) {
      sigmay = 0.2;
    }
    
    printf ("Gaussiana:\n\tx) N(%g, %g), y) N(%g, %g)\n\n", 
	    mediax, sigmax, mediay, sigmay);
    
    gx = newGaussiana (20, mediax, sigmax, 0);
    gy = newGaussiana (20, mediay, sigmay, 0);
    
    imprimirLog ("reconstructor.log", "Parametros libres iniciales: \n");
    for (i = 0; i < r->fL->n_pols; i++) {
      r->p[3 * i]     = gauss (gx);
      r->p[3 * i + 1] = gauss (gy);
      r->p[3 * i + 2] = init_value;
      imprimirLog ("reconstructor.log", "%g\t%g\t%g\n", 
		   r->p[3 * i], r->p[3 * i + 1], r->p[3 * i + 2]);
    }
    
    free (gx);
    free (gy);
  }
  else {    
    imprimirLog ("reconstructor.log", "inicializando con distribucion uniforme\n");
    imprimirLog ("reconstructor.log", "Parametros libres iniciales: \n");
    for (i = 0; i < r->fL->n_pols; i++) {
      r->p[3 * i]     = (float) random() / RAND_MAX;
      r->p[3 * i + 1] = (float) random() / RAND_MAX;
      r->p[3 * i + 2] = init_value;
      imprimirLog ("reconstructor.log", "%g\t%g\t%g\n", 
		   r->p[3 * i], r->p[3 * i + 1], r->p[3 * i + 2]);
    }
  }
  return r;
}

double run (Reconstructor *r, double ftol) {

  int i, j, n, nVis, status, seed = 0, npikaia;
  double func = 1e299, func_old = 1e300, x, y, valor, 
    funcMin = 1e300, init_value = r->p[2];
  float *p_old = NULL, *p_min;
  FILE *archivo = fopen("L_n.dat", "w");
  double *parsD, cuantaSize = 0;
  float ctrl[12];

  fclose(archivo);
  printf ("oki\n");
  nVis = r->nVis;
  printf ("okii\n");

  n = r->fL->n_pols;
  r->ftol = ftol;

  // itera sobre el numero de centros de Voronoi de la imagen.
  
  parsD = (double *) malloc (3 * n * sizeof (double));
  for (i = 0; i < r->fL->n_pols; i++) {
    parsD[3 * i]     = (double) r->p[3 * i];
    parsD[3 * i + 1] = (double) r->p[3 * i + 1];
    parsD[3 * i + 2] = (double) r->p[3 * i + 2];
  }

  func = L (r->fL, parsD, n);
  free (parsD);
  do_write_fits(r->fL->fg_image, "!func1.fits");
  printf ("func = %g, n = %d\n", func, n);
  imprimirLog ("reconstructor.log", "run:\n");

  imprimirLog ("reconstructor.log", "    parametros libres:\n");
  for (i = 0; i < n; i++) {
    printf ("  p[%d] = (%g, %g, %g)\n", i, r->p[3 * i], 
	    r->p[3 * i + 1], r->p[3 * i + 2]);
    imprimirLog ("reconstructor.log", " p[%d] = (%g, %g, %g)\n", i, r->p[3 * i], 
		 r->p[3 * i + 1], r->p[3 * i + 2]);
  }

  imprimirLog ("reconstructor.log", "    n = %d\n", n);
  archivo = fopen ("L_n.dat", "a");
  func_old = func;
  
  
  if (r->fL != NULL) {
    cuantaSize = r->fL->difmapNoise;
    eliminarFuncL (r->fL);
    r->fL = NULL;
  }
  r->fL = newFuncL (r->nombreImagen, r->nombresVis, nVis, n, 
		    r->entropia, cuantaSize);
  
  r->iter = 0;

  ctrl [0] = r->individuos;   // individuos (100)
  ctrl [1] = r->generaciones; // generaciones (500)
  ctrl [2] = -1;  // digitos significativos (6)
  ctrl [3] = -1;  // pbb de cruza (0.85)
  ctrl [4] = -1;  // modo de mutacion (2)
  ctrl [5] = -1;  // rango de mutacion inicial (0.005)
  ctrl [6] = -1;  // rango de mutacion minimo (0.0005)
  ctrl [7] = -1;  // rango de mutacion maximo (0.25)
  ctrl [8] = -1;  // diferencial relativo de ajuste (1)
  ctrl [9] = -1;  // plan de reproduccion (3)
  ctrl [10] = -1; // bandera de elitismo (0)
  ctrl [11] = -1; // imprimir output (0)

  rninit_(&seed);
  
  r_out = r;
  npikaia = 3 * n;
  pikaia_(ff_, &npikaia, ctrl, r->p, &func, &status);
  //func =  minimizador (r, ftol);

  imprimirLog ("reconstructor.log", "      func = %g\n", func);
  //func = L (r->fL, r->p, n + 1);
  printf ("funcOld = %g, func = %g, x = %g, y= %g, valor = %g, n = %d\n",
	  func_old, func, x, y, valor, n);

  
  if (r->entropia) {
    fprintf(archivo, "%d\t%lf\n", n, func);
  }
  else {
    int arch, nVis = 0;
    for (arch = 0; arch < r->fL->n_archivos; arch++) {
      nVis += r->fL->header_obs[arch]->nif * r->fL->header_obs[arch]->nsamp;
    }
    fprintf(archivo, "%d\t%lf\t%lf\n", n, func, 2 * func/ (2 * nVis - 3 * n));
    
    
    //printf("2) x = %g, y = %g, valor = %g\n", x, y, valor);
    //}
  
    fclose (archivo);
  }
  
  printf ("FIN, func_old = %.20g, func = %.20g\n", func_old, func);
  archivo = fopen ("L_n.dat", "a");
  fprintf (archivo, "FIN, func_old = %.20g, func = %.20g\n", func_old, func);
  fclose (archivo);
  free (r->p);
  r->p = p_old;
  return func;
}

void calcularNuevaPosicion(Reconstructor *r, double *x, double *y, 
			   double *valor) {
  int i, j, id, *ids, n = 1e4, ni = 0;
  double Itotal = 0, Imin = 0;
  
  
  for (i = 0; i < r->fL->n_pols; i++) {
    if (r->p[3 * i + 2] < Imin) {
      Imin = r->p[3 * i + 2];
    }
  }
  
  for (i = 0; i < r->fL->n_pols; i++) {
    Itotal += r->p[3 * i + 2] - Imin;
  }
  ids = (int *) malloc (n * sizeof (int));
  for (i = 0; i < r->fL->n_pols; i++) {
    for (j = 0; j < round (n * (r->p[3 * i + 2] - Imin) / Itotal); j++) {
      if (ni >= n) {
	break;
      }
      ids[ni] = i + 3;
      ni++;
    }
  }

  if (ni > n) {
    fprintf (stderr, "ERROR en calcularNuevaPosicion, (ni = %d) != (n = %d)\n", ni, n);
    exit (1);
  }

  id = ids[(int) floor(n * (float) random() / RAND_MAX)];
  free (ids);
  ids = NULL;
  if (id > r->fL->n_pols + 3) {
    fprintf (stderr, "ERROR en calcularNuevaPosicion, id > %d\n", r->fL->n_pols);
    exit (1);
  }
  (*valor) = r->p[3 * (id - 3) + 2];
  
  n = 0;
  for (i = 0; i < r->fL->fg_image->size[0]; i++) {
    for (j = 0; j < r->fL->fg_image->size[1]; j++) {
      if (r->fL->mask[i][j] == id) {
	n++;
      }
    }
  }
  
  n = (int) floor(n * (float) random() / RAND_MAX);
  ni = 0;
  for (i = 0; i < r->fL->fg_image->size[0]; i++) {
    for (j = 0; j < r->fL->fg_image->size[1]; j++) {
      if (r->fL->mask[i][j] == id) {
	ni++;
	if (ni == n) {
	  (*x) = (double) i / (r->fL->fg_image->size[0] - 1);
	  (*y) = (double) j / (r->fL->fg_image->size[1] - 1);
	  /*	  if (!validarMallaPoligonosAristas (r->fL->malla)) {
	    fprintf (stderr, "ERROR al validar malla en CalcularNuevaPosicion\n");
	    exit (1);
	    }*/
	  return;
	}
      }
    }
  }

  (*x) = 0.5;
  (*y) = 0.5;
}

double calcularPromedio (MallaVoronoi *m, double x, double y) {
  PoligonoVoronoi *pol;
  AristaVoronoi *a, *a_old;
  double suma = 0;
  int cont = 2;

  pol = encontrarPoligono (m, x, y);

  suma = pol->valor;

  a_old = pol->a;
  if (a_old->poliDer == pol) {
    suma += a_old->poliIzq->valor;
    a = a_old->ccwSucc;
  }
  else if (a_old->poliIzq == pol) {
    suma += a_old->poliIzq->valor;
    a = a_old->ccwPred;
  }
  else {
    fprintf (stderr, "ERROR al calcular promedio\n.");
    exit (1);
    //return 0;
  }

  while (a != a_old) {
    cont++;
    if (a->poliDer == pol) {
      suma += a->poliIzq->valor;
      a = a->ccwSucc;
    }
    else if (a->poliIzq == pol) {
      suma += a->poliDer->valor;
      a = a->ccwPred;
    }
    else {
      fprintf (stderr, "ERROR al calcular promedio\n.");
      exit (1);
    }
    
  }
  
  return suma / cont;
}

/*********************************************************************
 * Cambia la posicion de los poligonos que caigan en una zona donde la
 * atenuacion sea menor que r->Acutoff.
 ********************************************************************/

void truncarPorAtenuacion (Reconstructor *r){
  int i, j, k, nx, ny;
  float x, y;

  nx = r->fL->fg_image->size[0];
  ny = r->fL->fg_image->size[1];

  for (k = 0; k < r->fL->n_pols; k++) {
    x = r->p [k * 3];	 
    y = r->p [k * 3 + 1];
    i = (int) round (x * (nx - 1.0));
    j = (int) round (y * (ny - 1.0));
    if (i < 0 || i >= ny || j < 0 || j >= ny || 
	r->fL->atten[0][i + j * r->fL->fg_image->size[0] + 1][0] < r->Acutoff) {
      printf (">>>>>>>>>Encontrando para (%d, %d) = (%g, %g)\n", i, j, 
	      x, y);
      encontrarPosicionAtenuacion (r, &i, &j, 1);
      printf ("Cambiando (%f, %f) por (%g, %g)\n", 
	      r->p [k * 3], r->p [k * 3 + 1], i / (nx - 1.0), j / (ny - 1.0));
      r->p [k * 3]     = i / (nx - 1.0);
      r->p [k * 3 + 1] = j / (ny - 1.0);
      r->p [k * 3 + 2] = r->fL->fg_image->pixels[i + j * nx] / r->fL->difmapNoise;
    }
  }
}

/*********************************************************************
 * Busca el punto mas cercano a (i_ori, j_ori) cuya atenuacion sea
 * mayor a r->Acutoff.
 ********************************************************************/

void encontrarPosicionAtenuacion (Reconstructor *r, int *i_ori, int *j_ori,
				  int nivel) {
  float dist, dist_min = 1e10;
  int i, j, i_final, j_final, nx, ny;

  nx = r->fL->fg_image->size[0];
  ny = r->fL->fg_image->size[1];

  j = (*j_ori) + nivel;
  if (j >= 0 && j < ny) {
    for (i = (*i_ori) - nivel; i <= (*i_ori) + nivel; i++) {
      if (i >= 0 && i < nx &&
	  r->fL->atten[0][i + j * r->fL->fg_image->size[0] + 1][0] >= r->Acutoff) {
	dist = (i - (*i_ori))*(i - (*i_ori)) + (j - (*j_ori))*(j - (*j_ori));
	if (dist < dist_min) {
	  i_final = i;
	  j_final = j;
	  dist_min = dist;
	}
      }
    }
  }
  j = (*j_ori) - nivel;
  if (j >= 0 && j < ny) {
    for (i = (*i_ori) - nivel; i <= (*i_ori) + nivel; i++) {
      if (i >= 0 && i < nx &&
	  r->fL->atten[0][i + j * r->fL->fg_image->size[0] + 1][0] >= r->Acutoff) {
	dist = (i - (*i_ori))*(i - (*i_ori)) + (j - (*j_ori))*(j - (*j_ori));
	if (dist < dist_min) {
	  i_final = i;
	  j_final = j;
	  dist_min = dist;
	}
      }
    }
  }

  i = (*i_ori) + nivel;
  if (i >= 0 && i < nx) {
    for (j = (*j_ori) - nivel; j <= (*j_ori) + nivel; j++) {
      if (j >= 0 && j < ny &&
	  r->fL->atten[0][i + j * r->fL->fg_image->size[0] + 1][0] >= r->Acutoff) {
	dist = (i - (*i_ori))*(i - (*i_ori)) + (j - (*j_ori))*(j - (*j_ori));
	if (dist < dist_min) {
	  i_final = i;
	  j_final = j;
	  dist_min = dist;
	}
      }
    }
  }
  i = (*i_ori) + nivel;
  if (i >= 0 && i < nx) {
    for (j = (*j_ori) - nivel; j <= (*j_ori) + nivel; j++) {
      if (j >= 0 && j < ny &&
	  r->fL->atten[0][i + j * r->fL->fg_image->size[0] + 1][0] >= r->Acutoff) {
	dist = (i - (*i_ori))*(i - (*i_ori)) + (j - (*j_ori))*(j - (*j_ori));
	if (dist < dist_min) {
	  i_final = i;
	  j_final = j;
	  dist_min = dist;
	}
      }
    }
  }

  if (dist_min == 1e10) {
    encontrarPosicionAtenuacion (r, i_ori, j_ori, nivel + 1);
  }
  else {
    (*i_ori) = i_final;
    (*j_ori) = j_final;
  }
}

void reinicializar (Reconstructor *r, char *nombreArchivo){
  FILE *archivo = fopen (nombreArchivo, "r");
  int n, i;
  double x, y, valor, cuantaSize = 0;

  printf ("reinicializando con %s\n", nombreArchivo);
  fscanf (archivo, "%d\n", &n);
  if (r->fL != NULL) {
    cuantaSize = r->fL->difmapNoise;
    eliminarFuncL (r->fL);
    r->fL = NULL;
  }
  /*  else {
    eliminarMallaVoronoi (r->fL->malla);
    r->fL->n_pols = n;
    r->fL->malla = newMallaVoronoi();
    }*/
  r->fL = newFuncL (r->nombreImagen, r->nombresVis, r->nVis, n, 
		    r->entropia, r->fL->difmapNoise);
  r->fL->malla = newMallaVoronoi ();
  if (r->p != NULL) {
    free (r->p);
    r->p = NULL;
  }
  r->p = (float *) malloc (3 * r->fL->n_pols * sizeof (float));
  for (i = 0; i < n; i++) {
    fscanf (archivo, "%lf\t%lf\t%lf\n", &x, &y, &valor);
    insertarSitio (r->fL->malla, x, y, valor);
    r->p[3 * i]     = x;
    r->p[3 * i + 1] = y;
    r->p[3 * i + 2] = valor / r->fL->difmapNoise;
  }
  fclose (archivo);
  printf ("reinicializado con %s\n", nombreArchivo);
}

void eliminarReconstructor (Reconstructor *r) {
  int i;

  free (r->nombreImagen);
  for (i = 0; i < r->nVis; i++){
    free (r->nombresVis[i]);
  }
  free (r->nombresVis);
  eliminarFuncL (r->fL);
  free (r);
}

/*********************************************************
   Retorna un nuevo reconstructor leyendo los parametros 
   desde el arhivo de entrada nombreArchivo.
*********************************************************/

Reconstructor *leerArchivoEntrada (char *nombreArchivo) {
  FILE *archivo = fopen (nombreArchivo, "r");
  char *nombreImagen, **nombresVisibilidades, *dato;
  int nPols, entropia, nVis, i, init_gauss;
  double init_value, cuantaSize, cutoff, Nmax;
  Reconstructor *r;
  float individuos, generaciones;

  if (archivo == NULL) {
    fprintf (stderr, "ERROR al intentar leer %s\n", nombreArchivo);
    exit (1);
  }
  
  nombreImagen = (char *) malloc (512 * sizeof (char));
  dato = (char *) malloc (512 * sizeof (char));
  
  fscanf (archivo, "%s\t%s\n", dato, nombreImagen);
  printf ("%s\t%s\n", dato, nombreImagen);
  fscanf (archivo, "%s\t%d\n", dato, &nPols);
  printf ("%s\t%d\n", dato, nPols);
  fscanf (archivo, "%s\t%d\n", dato, &entropia);
  printf ("%s\t%d\n", dato, entropia);
  fscanf (archivo, "%s\t%lf\n", dato, &init_value);
  printf ("%s\t%g\n", dato, init_value);
  fscanf (archivo, "%s\t%d\n", dato, &init_gauss);
  printf ("%s\t%d\n", dato, init_gauss);
  fscanf (archivo, "%s\t%lf\n", dato, &cuantaSize);
  printf ("%s\t%g\n", dato, cuantaSize);
  fscanf (archivo, "%s\t%lf\n", dato, &cutoff);
  printf ("%s\t%g\n", dato, cutoff);
  fscanf (archivo, "%s\t%lf\n", dato, &Nmax);
  printf ("%s\t%g\n", dato, Nmax);
  fscanf (archivo, "%s\t%f\n", dato, &individuos);
  printf ("%s\t%g\n", dato, individuos);
  fscanf (archivo, "%s\t%f\n", dato, &generaciones);
  printf ("%s\t%g\n", dato, generaciones);
  fscanf (archivo, "%s\t%d\n", dato, &nVis);
  printf ("%s\t%d\n", dato, nVis);
 
  nombresVisibilidades = (char **) malloc (nVis * sizeof(char *));
  for (i = 0; i < nVis; i++) {
    nombresVisibilidades[i] = (char *) malloc (512 * sizeof (char));
    fscanf (archivo, "%s\n", nombresVisibilidades[i]);
    printf ("%s\n", nombresVisibilidades[i]);
  }

  r = newReconstructor (nombreImagen, nombresVisibilidades, 
			nVis, nPols, init_value, init_gauss, entropia, 
			cuantaSize, cutoff);
  r->Nmax = Nmax;
  r->individuos = individuos;
  r->generaciones = generaciones;
  fclose (archivo);
  free (nombreImagen);
  free (dato);
  for (i = 0; i < nVis; i++) {
    free (nombresVisibilidades[i]);
  }
  printf ("leer: r->nVis = %d\n\n", r->nVis);

  return r;
}

void imprimirLog (char *archivo, char* s, ...) {
  va_list ap;
  FILE* archivoLog = fopen (archivo, "a");

  va_start (ap, s);
  vfprintf (archivoLog, s, ap);
  va_end (ap);

  fclose (archivoLog);
}

/*
void imprimirPikaia (char* s, ...) {
  va_list ap;
  FILE* archivoLog = fopen ("pikaia.log", "a");

  va_start (ap, s);
  vfprintf (archivoLog, s, ap);
  va_end (ap);

  fclose (archivoLog);
  }*/
