#include "reconstructor_aristas.h"

Reconstructor *newReconstructor (char * nombreImagen, char **nombresVis, 
				 int nVis, int n,
				 double init_value, int init_gauss,
				 int entropia, double cuantaSize) {
  int i, j, k, nx, ny, imax = 0, jmax = 0;
  FILE* archivoLog = fopen ("reconstructor.log", "w");
  fclose (archivoLog);

  Reconstructor *r = (Reconstructor *) malloc (sizeof (Reconstructor));

  r->entropia = entropia;
  r->nVis = nVis;
  printf ("nVis = %d\n", nVis);
  printf ("r->nVis = %d\n\n", r->nVis);
  r->iter = 0;
  r->fL = newFuncL(nombreImagen, nombresVis, nVis, n, r->entropia, cuantaSize);
  r->nombreImagen = (char *) malloc (256 * sizeof (char));
  strcpy(r->nombreImagen, nombreImagen);
  r->nombresVis = (char **) malloc (nVis * sizeof (char *));
  imprimirLog ("Reconstruyendo ");
  for (i = 0; i < nVis; i++){
    r->nombresVis[i] = (char *) malloc (256 * sizeof (char));
    strcpy(r->nombresVis[i], nombresVis[i]);
    imprimirLog ("%s ", r->nombresVis[i]);
  }
  // Inicializamos la Imagen

  r->p = (float *) malloc (3 * r->fL->n_pols * sizeof (float));
  r->iter = 0;

  nx = r->fL->fg_image->size[0];
  ny = r->fL->fg_image->size[1];

  imprimirLog ("\nutilizando como imagen de referencia%s\n ",
	       r->nombreImagen);
  imprimirLog ("nx = %d, ny = %d\n ", nx, ny);
  imprimirLog ("entropia: %d\n", r->entropia);
  if (init_gauss) {
    double mediax, sigmax, mediay, sigmay, atten_max = -1e300;
    Gaussiana *gx, *gy;
    imprimirLog ("inicializando con distribucion gaussiana\n");

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
    
    imprimirLog ("Parametros libres iniciales: \n");
    for (i = 0; i < r->fL->n_pols; i++) {
      r->p[3 * i]     = gauss (gx);
      r->p[3 * i + 1] = gauss (gy);
      r->p[3 * i + 2] = init_value;
      imprimirLog ("%g\t%g\t%g\n", 
		   r->p[3 * i], r->p[3 * i + 1], r->p[3 * i + 2]);
    }
    
    free (gx);
    free (gy);
  }
  else {    
    imprimirLog ("inicializando con distribucion uniforme\n");
    imprimirLog ("Parametros libres iniciales: \n");
    for (i = 0; i < r->fL->n_pols; i++) {
      r->p[3 * i]     = (float) random() / RAND_MAX;
      r->p[3 * i + 1] = (float) random() / RAND_MAX;
      r->p[3 * i + 2] = init_value;
      imprimirLog ("%g\t%g\t%g\n", 
		   r->p[3 * i], r->p[3 * i + 1], r->p[3 * i + 2]);
    }
  }
  return r;
}

/*********************************************************
   Minimiza utilizando Polak-Ribiere de Numerical recipes.
   Retorna el valor final de la funcion a minimizar.
*********************************************************/

double minimizador (Reconstructor *r, double ftol) {

  float f (float pars[]){
    double ret; 
    double *parsD = (double *) malloc (3 * r->fL->n_pols * sizeof (double));
    int i, j;
    
    imprimirLog ("        Calculando funcion con parametros:\n");
    for (i = 0; i < r->fL->n_pols; i++) {
      parsD[3 * i]     = (double) pars[3 * i + 1];
      parsD[3 * i + 1] = (double) pars[3 * i + 1 + 1];
      // truncamos las intensidades.
      if (pars[3 * i + 2 + 1] <= 0) {
	pars[3 * i + 2 + 1] = MIN_PIX;
      }
      parsD[3 * i + 2] = (double) pars[3 * i + 2 + 1];
      imprimirLog ("        %g\t%g\t%g\n", 
		   parsD[3 * i], parsD[3 * i + 1], parsD[3 * i + 2]);
    }
    
    //printf ("Calculando chi2\n");
    ret  =  L(r->fL, parsD, r->fL->n_pols);
    printf ("chi2 = %g\n", ret);
    imprimirLog ("        func = %g\n\n", ret);
    free (parsD);
    /*
    if (r->iter <= 1) {
      mask = do_read (r->fL->nombreFits);
      for (i = 0; i < mask->size[0]; i++) {
	for (j = 0; j < mask->size[0]; j++) {
	  mask->pixels[i + j * mask->size[0]] = r->fL->mask[i][j];
	  //printf ("mask[%d][%d] = %g\n", i, j, fL->mask[i][j]);
	}
      }
      do_write_fits(mask, "!mask.fits");
      }*/
    
    return ret;
  }
  
  void df (float pars[], float grad[]) {
    double *parsD = (double *) malloc (3 * r->fL->n_pols * sizeof (double));
    double *gradD = (double *) malloc (3 * r->fL->n_pols * sizeof (double));
    int i, j;
    char fileout[30];
    struct image *imagenPix;
    struct image *mask;
    FILE *archivo;
    
    r->iter ++;
    imprimirLog ("         Calculando gradiente con:\n");
    for (i = 0; i < r->fL->n_pols; i++) {
      parsD [3 * i]     = (double) pars[3 * i + 1];
      parsD [3 * i + 1] = (double) pars[3 * i + 1 + 1];
      // truncamos las intensidades.
      if (pars[3 * i + 2 + 1] <= 0) {
	pars[3 * i + 2 + 1] = MIN_PIX;
      }
      parsD [3 * i + 2] = (double) pars[3 * i + 2 + 1];
      imprimirLog ("        %g\t%g\t%g\n", 
		   parsD[3 * i], parsD[3 * i + 1], parsD[3 * i + 2]);
    }
    printf ("Calculando dchi2\n");
    dL(r->fL, parsD, gradD, r->fL->n_pols, 0);
    
    //imprimirLog ("\n");
    sprintf(fileout, "gradiente_%d_%d.fits", r->fL->n_pols, r->iter);
    archivo = fopen (fileout, "w");
    imprimirLog ("Gradiente:\n");
    for (i = 0; i < r->fL->n_pols; i++) {
      grad [3 * i + 1] = (float) gradD [3 * i];
      grad [3 * i + 1 + 1] = (float) gradD [3 * i + 1];
      grad [3 * i + 2 + 1] = (float) gradD [3 * i + 2];
      imprimirLog ("        grad[%d] = %g, grad[%d] = %g, grad[%d] = %g\n",
		   3 * i + 1, grad [3 * i + 1], 
		   3 * i + 1 + 1, grad [3 * i + 1 + 1], 
		   3 * i + 1 + 2, grad [3 * i + 1 + 2]);
      printf ("grad[%d] = %g, grad[%d] = %g, grad[%d] = %g\n",
	      3 * i + 1, grad [3 * i + 1], 
	      3 * i + 1 + 1, grad [3 * i + 1 + 1], 
	      3 * i + 1 + 2, grad [3 * i + 1 + 2]);
      fprintf (archivo, "%d\t%g\t%g\t%g\n", i,
	       grad [3 * i + 1], grad [3 * i + 1 + 1], 
	       grad [3 * i + 1 + 2]);
    }
    imprimirLog ("\n");
    fclose(archivo);

    sprintf(fileout, "!mask_%d_%d.fits", r->fL->n_pols, r->iter);
    mask = do_read (r->fL->nombreFits);
    for (i = 0; i < mask->size[0]; i++) {
      for (j = 0; j < mask->size[0]; j++) {
	mask->pixels[i + j * mask->size[0]] = r->fL->mask[i][j];
	//printf ("mask[%d][%d] = %g\n", i, j, fL->mask[i][j]);
      }
    }
    do_write_fits(mask, fileout);

    sprintf(fileout, "!MEM_%d_%d.fits", r->fL->n_pols, r->iter);
    do_write_fits(r->fL->fg_image, fileout);
    
    // Imprimimos los pixeles integrados en la 1era iteracion. -------
    if (r->iter == 1) {
      int suma = 0;
      imagenPix = do_read(r->fL->nombreFits);
      for (i = 0; i < r->fL->fg_image->size[0]; i++) {
	for (j = 0; j < r->fL->fg_image->size[1]; j++) {
	  imagenPix->pixels[i + j * r->fL->fg_image->size[0]] = 
	    r->fL->PixIntegrados_x[i][j];
	  suma += r->fL->PixIntegrados_x[i][j];
	}
      }
      sprintf (fileout, "!PixIntegrados_x.fits");
      do_write_fits (imagenPix, fileout);
      delete_map (imagenPix);
      printf ("suma = %d\n\n", suma);
      
      imagenPix = do_read(r->fL->nombreFits);
      for (i = 0; i < r->fL->fg_image->size[0]; i++) {
	for (j = 0; j < r->fL->fg_image->size[1]; j++) {
	  imagenPix->pixels[i + j * r->fL->fg_image->size[0]] = 
	    r->fL->PixIntegrados_I[i][j];
	}
      }
      sprintf (fileout, "!PixIntegrados_I.fits");
      do_write_fits(imagenPix, fileout);
      delete_map (imagenPix);    
    }
    // ------------------------------------------------------
    
    // Imprimimos en un archivo la malla.
    sprintf(fileout, "malla_%d_%d.dat", r->fL->n_pols, r->iter);
    printf ("Imprimiendo %s\n", fileout);
    //imprimirMallaArchivo (r->fL->malla, fileout);
    printf ("Impresa\n");
    
    free (parsD);
    free (gradD);
    
    printf ("dchi2 calculado\n");
  }

  float fret;
  long t;
  FILE *archivo;
  char fileout[30];

  if (r->fL->n_pols < 10) {
    sprintf(fileout, "chi2_0%d.dat", r->fL->n_pols);
  }
  else {
    sprintf(fileout, "chi2_%d.dat", r->fL->n_pols);
  }
  archivo = fopen (fileout, "w");
 
  t = time(0);
  frprmn(r->p - 1, 3 * r->fL->n_pols, ftol, &(r->iter), &fret, f, df);
  t = time(0) - t;
  printf("Terminado P-R en %d segundos y %d iteraciones\n", 
	 (int) t, r->iter);

  if (r->fL->n_pols < 10) {
    sprintf(fileout, "!MEM_NR_00%d.fits", r->fL->n_pols);
  }
  else if (r->fL->n_pols < 100) {
    sprintf(fileout, "!MEM_NR_0%d.fits", r->fL->n_pols);
  }
  else {
    sprintf(fileout, "!MEM_NR_%d.fits", r->fL->n_pols);
  }
  do_write_fits(r->fL->fg_image, fileout);
  fprintf(archivo,
	  "Terminado P-R en %d segundos y %d iteraciones\n", 
	  (int) t, r->iter);  
  fprintf(archivo, "Chi2 = %.25g\n", fret);

  if (r->fL->n_pols < 10) {
    sprintf(fileout, "MEM_NR_00%d.dat", r->fL->n_pols);
  }
  else if (r->fL->n_pols < 100) {
    sprintf(fileout, "MEM_NR_0%d.dat", r->fL->n_pols);
  }
  else {
    sprintf(fileout, "MEM_NR_%d.dat", r->fL->n_pols);
  }
  imprimirMallaArchivo (r->fL->malla, fileout);
  fclose (archivo);
  //free (p);

  return fret;
}

double run (Reconstructor *r, double ftol) {
  int i, j, n, nVis;
  double func = 1e299, func_old = 1e300, x, y, valor, 
    funcMin = 1e300, init_value = r->p[2], cuantaSize;
  float *p_old = NULL, *p_min;
  FILE *archivo = fopen("L_n.dat", "w");
  double *parsD;

  fclose(archivo);
  printf ("oki\n");
  nVis = r->nVis;
  printf ("okii\n");

  n = r->fL->n_pols;
  /*
  // Buscamos una buena aproximacion para los poligonos iniciales.
  p_min = (float *) malloc (3 * n * sizeof (float));
  imprimirLog ("run: Buscando buena aproximacion inicial para n = %d.\n", n);
  for (j = 0; j < 100; j++) {
    printf ("---Buscando buena aproximacion inicial %d.\n", j);
    if (r->fL != NULL) {
      eliminarFuncL (r->fL);
      r->fL = NULL;
    }
    r->fL = newFuncL(r->nombreImagen, r->nombresVis, nVis, n, r->entropia,
		     r->fL->difmapNoise);
    imprimirLog ("  Parametros:\n");
    for (i = 0; i < n; i++) {
      imprimirLog ("    %g\t%g\t%g\n", r->p[3 * i], r->p[3 * i + 1], r->p[3 * i + 2]);
    }
    func =  minimizador (r, ftol);
    printf ("---func = %g\n", func);
    imprimirLog ("  func = %g\n\n", func);
    if (func < funcMin) {
      funcMin = func;
      for (i = 0; i < n; i++) {
	p_min[3 * i]     = r->p[3 * i];
	p_min[3 * i + 1] = r->p[3 * i + 1];
	p_min[3 * i + 2] = r->p[3 * i + 2];
      }
    }
    for (i = 0; i < n; i++) {
      r->p[3 * i]     = (float) random() / RAND_MAX;
      r->p[3 * i + 1] = (float) random() / RAND_MAX;
      r->p[3 * i + 2] = init_value;
    }
    printf ("---funcMin = %g\n", funcMin);
    imprimirLog ("  funcMin = %g\n", funcMin);
  }
  for (i = 0; i < n; i++) {
    r->p[3 * i]     = p_min[3 * i];
    r->p[3 * i + 1] = p_min[3 * i + 1];
    r->p[3 * i + 2] = p_min[3 * i + 2];
  }
  free (p_min);
  */

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
  imprimirLog ("run:  Primera evaluacion\n");

  for (i = 0; i < n; i++) {
    printf ("  p[%d] = (%g, %g, %g)\n", i, r->p[3 * i], 
	    r->p[3 * i + 1], r->p[3 * i + 2]);
    imprimirLog (" p[%d] = (%g, %g, %g)\n", i, r->p[3 * i], 
		 r->p[3 * i + 1], r->p[3 * i + 2]);
  }
  imprimirLog ("  func = %g\n", func);
  //for (n = r->fL->n_pols; func < func_old; n++) {
  imprimirLog ("  Comienzan iteraciones sobre el numero de poligonos\n");
  for (n = r->fL->n_pols + 1; n < 3000; n++) {
    MallaVoronoi *malla_old = newMallaVoronoi();

    imprimirLog ("    n = %d\n", n);
    archivo = fopen ("L_n.dat", "a");
    func_old = func;

    if (p_old != NULL) {
      free (p_old);
      p_old = NULL;
    }
    p_old = r->p;

    r->p = (float *) malloc (3 * n * sizeof (float));
    imprimirLog ("    parametros libres:\n");
    for (i = 0; i < n - 1; i++) {
      r->p[3 * i]     = p_old[3 * i];
      r->p[3 * i + 1] = p_old[3 * i + 1];
      r->p[3 * i + 2] = p_old[3 * i + 2];

      if (r->p[3 * i] < 0) {
	r->p[3 * i] = 0;
      }
      else if (r->p[3 * i] > 1) {
	r->p[3 * i] = 1;
      }
      if (r->p[3 * i + 1] < 0) {
	r->p[3 * i + 1] = 0;
      }
      else if (r->p[3 * i + 1] > 1) {
	r->p[3 * i + 1] = 1;
      }
      imprimirLog ("      %g\t%g\t%g\n", 
		   r->p[3 * i], r->p[3 * i + 1], r->p[3 * i + 2]);
      insertarSitio (malla_old, r->p[3 * i], r->p[3 * i + 1], 
		     r->p[3 * i + 2] * r->fL->difmapNoise);
    }
    imprimirLog ("\n");        

    calcularNuevaPosicionVertice (r, &x, &y, &valor);
    cuantaSize = r->fL->difmapNoise;
    if (r->fL != NULL) {
      eliminarFuncL (r->fL);
      r->fL = NULL;
    }
    r->fL = newFuncL(r->nombreImagen, r->nombresVis, nVis, n, r->entropia, 
		     cuantaSize);
    //do{
      //valor = calcularPromedio (r, x, y);
      //x = (float) random() / RAND_MAX;;
      //y = (float) random() / RAND_MAX;;
      //valor = calcularPromedio (malla_old, x, y) / r->fL->difmapNoise;
      r->p[3 * (n - 1)]     = x;
      r->p[3 * (n - 1) + 1] = y;
      r->p[3 * (n - 1) + 2] = valor;
      //for (i = 0; i < n; i++) {
      //r->p[3 * i + 2] = 0;
      //}
      imprimirLog ("      posible nueva posicion: %g\t%g\t%g\n", x, y, valor);

      r->iter = 0;
      func =  minimizador (r, ftol);
      imprimirLog ("      func = %g\n", func);
      //func = L (r->fL, r->p, n + 1);
      printf ("funcOld = %g, func = %g, x = %g, y= %g, valor = %g, n = %d\n",
	      func_old, func, x, y, valor, n);
      //}while (func > func_old);
    imprimirLog ("    nueva posicion definitiva: %g\t%g\t%g\n", x, y, valor);

    if (r->entropia) {
      fprintf(archivo, "%d\t%lf\n", n, func);
    }
    else {
      int arch, nVis = 0;
      for (arch = 0; arch < r->fL->n_archivos; arch++) {
	nVis += r->fL->header_obs[arch]->nif * r->fL->header_obs[arch]->nsamp;
      }
      fprintf(archivo, "%d\t%lf\t%lf\n", n, func, 2 * func/ (2 * nVis - 3 * n));
    }

    
      //printf("2) x = %g, y = %g, valor = %g\n", x, y, valor);
      //}
      
    fclose (archivo);
    eliminarMallaVoronoi (malla_old);
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
  //int i, j, id, *ids, n = 1e4, ni = 0;
  double dITotal = 0, I;
  NodoLista *n;
  AristaVoronoi *a;
  
  for (n = r->fL->malla->aristas->primero; n != NULL; n = n->sgte) {
    a = (AristaVoronoi *) n->info;
    //imprimirArista (a);
    if (a->poliDer != NULL && a->poliIzq != NULL &&
	a->poliDer->id >= 3 && a->poliIzq->id >= 3) {
      dITotal += fabs (a->poliDer->valor - a->poliIzq->valor);
    }
  }

  I = dITotal * random () / RAND_MAX;

  dITotal = 0;
  for (n = r->fL->malla->aristas->primero; n != NULL; n = n->sgte) {
    a = (AristaVoronoi *) n->info;
    if (a->poliDer != NULL && a->poliIzq != NULL &&
	a->poliDer->id >= 3 && a->poliIzq->id >= 3) {
      dITotal += fabs (a->poliDer->valor - a->poliIzq->valor);
      
      if (dITotal > I) {
	(*x)     = (a->ptoIni->x + a->ptoFin->x) / 2;
	(*y)     = (a->ptoIni->y + a->ptoFin->y) / 2;
	(*valor) = (a->poliDer->valor + a->poliIzq->valor) / 2;
	return;
      }
    }
  }
  (*x)     = (a->ptoIni->x + a->ptoFin->x) / 2;
  (*y)     = (a->ptoIni->y + a->ptoFin->y) / 2;
  (*valor) = (a->poliDer->valor + a->poliIzq->valor) / 2;
  printf ("calcularNuevaPosicion: Se llego al final.\n");
}

void calcularNuevaPosicionVertice (Reconstructor *r, double *x, double *y, 
				   double *valor) {
  //int i, j, id, *ids, n = 1e4, ni = 0;
  double dITotal = 0, I, v;
  NodoLista *n;
  PuntoVoronoi *p;
  AristaVoronoi *a;
  
  for (n = r->fL->malla->puntos->primero; n != NULL; n = n->sgte) {
    p = (PuntoVoronoi *) n->info;
    a = p->a;
    //imprimirArista (a);
    if (a->poliDer != NULL && a->poliIzq != NULL && 
	a->poliDer->id >= 3 && a->poliIzq->id >= 3) {
      if (a->ptoIni == p){
	if (a->cwPred->ptoIni == p) {
	  if (a->cwPred->poliDer != NULL && a->cwPred->poliDer->id >= 3) {
	    dITotal += (fabs (a->poliDer->valor - a->poliIzq->valor) +
			fabs (a->cwPred->poliDer->valor - a->poliIzq->valor) +
			fabs (a->poliDer->valor - a->cwPred->poliDer->valor));
	  }
	}
	else {
	  if (a->cwPred->poliIzq != NULL && a->cwPred->poliIzq->id >= 3) {
	    dITotal += (fabs (a->poliDer->valor - a->poliIzq->valor) +
			fabs (a->cwPred->poliIzq->valor - a->poliIzq->valor) +
			fabs (a->poliDer->valor - a->cwPred->poliIzq->valor));
	  }
	}
      }
      else {
	if (a->ccwSucc->ptoIni == p){
	  if (a->ccwSucc->poliIzq != NULL && a->ccwSucc->poliIzq->id >= 3) {
	    dITotal += (fabs (a->poliDer->valor - a->poliIzq->valor) +
			fabs (a->ccwSucc->poliIzq->valor - a->poliIzq->valor) +
			fabs (a->poliDer->valor - a->ccwSucc->poliIzq->valor));
	  }
	}
	else {
	  if (a->ccwSucc->poliDer != NULL && a->ccwSucc->poliDer->id >= 3) {
	    dITotal += (fabs (a->poliDer->valor - a->poliIzq->valor) +
			fabs (a->ccwSucc->poliDer->valor - a->poliIzq->valor) +
			fabs (a->poliDer->valor - a->ccwSucc->poliDer->valor));
	  }
	}
      }
    }
  }

  I = dITotal * random () / RAND_MAX;

  dITotal = 0;
  for (n = r->fL->malla->puntos->primero; n != NULL; n = n->sgte) {
    p = (PuntoVoronoi *) n->info;
    a = p->a;
    //imprimirArista (a);
    if (a->poliDer != NULL && a->poliIzq != NULL && 
	a->poliDer->id >= 3 && a->poliIzq->id >= 3) {
      if (a->ptoIni == p){
	if (a->cwPred->ptoIni == p) {
	  if (a->cwPred->poliDer != NULL && a->cwPred->poliDer->id >= 3) {
	    dITotal += (fabs (a->poliDer->valor - a->poliIzq->valor) +
			fabs (a->cwPred->poliDer->valor - a->poliIzq->valor) +
			fabs (a->poliDer->valor - a->cwPred->poliDer->valor));
	    v = (a->poliDer->valor + a->poliIzq->valor + a->cwPred->poliDer->valor) / 3;
	  }
	}
	else {
	  if (a->cwPred->poliIzq != NULL && a->cwPred->poliIzq->id >= 3) {
	    dITotal += (fabs (a->poliDer->valor - a->poliIzq->valor) +
			fabs (a->cwPred->poliIzq->valor - a->poliIzq->valor) +
			fabs (a->poliDer->valor - a->cwPred->poliIzq->valor));
	    v = (a->poliDer->valor + a->poliIzq->valor + a->cwPred->poliIzq->valor) / 3;
	  }
	}
      }
      else {
	if (a->ccwSucc->ptoIni == p){
	  if (a->ccwSucc->poliIzq != NULL && a->ccwSucc->poliIzq->id >= 3) {
	    dITotal += (fabs (a->poliDer->valor - a->poliIzq->valor) +
			fabs (a->ccwSucc->poliIzq->valor - a->poliIzq->valor) +
			fabs (a->poliDer->valor - a->ccwSucc->poliIzq->valor));
	    v = (a->poliDer->valor + a->poliIzq->valor + a->ccwSucc->poliIzq->valor) / 3;
	  }
	}
	else {
	  if (a->ccwSucc->poliDer != NULL && a->ccwSucc->poliDer->id >= 3) {
	    dITotal += (fabs (a->poliDer->valor - a->poliIzq->valor) +
			fabs (a->ccwSucc->poliDer->valor - a->poliIzq->valor) +
			fabs (a->poliDer->valor - a->ccwSucc->poliDer->valor));
	    v = (a->poliDer->valor + a->poliIzq->valor + a->ccwPred->poliDer->valor) / 3;
	  }
	}
      }
    }
    if (dITotal > I) {
      (*x)     = p->x;
      (*y)     = p->y;
      (*valor) = v;
      return;
    }
  }
  (*x)     = p->x;
  (*y)     = p->y;
  (*valor) = v;
  printf ("calcularNuevaPosicion: Se llego al final.\n");
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

void reinicializar (Reconstructor *r, char *nombreArchivo){
  FILE *archivo = fopen (nombreArchivo, "r");
  int n, i;
  double x, y, valor;

  printf ("reinicializando con %s\n", nombreArchivo);
  fscanf (archivo, "%d\n", &n);
  if (r->fL != NULL) {
    eliminarFuncL (r->fL);
    r->fL = NULL;
  }
  /*  else {
    eliminarMallaVoronoi (r->fL->malla);
    r->fL->n_pols = n;
    r->fL->malla = newMallaVoronoi();
    }*/
  r->fL = newFuncL (r->nombreImagen, r->nombresVis, r->nVis, n, r->entropia,
		    r->fL->difmapNoise);
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
  eliminarFuncL(r->fL);
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
  double init_value, cuantaSize, cutoff;
  Reconstructor *r;

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
  fscanf (archivo, "%s\t%d\n", dato, &nVis);
  printf ("%s\t%d\n", dato, nVis);
 
  nombresVisibilidades = (char **) malloc (nVis * sizeof(char *));
  for (i = 0; i < nVis; i++) {
    nombresVisibilidades[i] = (char *) malloc (512 * sizeof (char));
    fscanf (archivo, "%s\n", nombresVisibilidades[i]);
    printf ("%s\n", nombresVisibilidades[i]);
  }

  r = newReconstructor (nombreImagen, nombresVisibilidades, 
			nVis, nPols, init_value, init_gauss, 
			entropia, cuantaSize);
  fclose (archivo);
  free (nombreImagen);
  free (dato);
  for (i = 0; i < nVis; i++) {
    free (nombresVisibilidades[i]);
  }
  printf ("leer: r->nVis = %d\n\n", r->nVis);

  return r;
}

void imprimirLog (char* s, ...) {
  va_list ap;
  FILE* archivoLog = fopen ("reconstructor.log", "a");

  va_start (ap, s);
  vfprintf (archivoLog, s, ap);
  va_end (ap);

  fclose (archivoLog);
}
