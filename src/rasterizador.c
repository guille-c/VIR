#include "rasterizador.h"

double **toImage(MallaVoronoi *malla, int nx, int ny, int **mask) {
  double **im = (double **) malloc(nx * sizeof(double *));
  int i, j;
  double m;
  NodoLista *n;
  AristaVoronoi *a;
  PoligonoVoronoi *polD, *polI;

  if (im == NULL) {
    fprintf (stderr, "ERROR en toImage, im = NULL.\n");
    exit (1);
  }
 
  //------------Rasterizador semi-bruto-------------------
  polD = (PoligonoVoronoi *) malla->poligonos->primero->info;
  for (i = 0; i < nx; i++) {
      im[i] = (double *) malloc(ny * sizeof(double));
    for (j = 0; j < ny; j++) {
      double x, y;
      x = i / (nx - 1.0);
      y = j / (ny - 1.0);
      
      polD = encontrarPoligono2 (polD, x, y);
      im [i][j] = polD->valor;
      if (mask != NULL) {
	mask[i][j] = polD->id;
      }
    }
  }
  
  return im;
  // -----------------------------------------------------

  if (mask == NULL) {
    for(i = 0; i < nx; i++) {
      im[i] = (double *) malloc(ny * sizeof(double));
      if (im[i] == NULL) {
	fprintf (stderr, "ERROR en toImage, im[%d] = NULL.\n", i);
	exit (1);
      }
      for (j = 0; j < ny; j++) {
	im[i][j] = NULL_PIX;
      }
    }
  }
  else {
    for(i = 0; i < nx; i++) {
      im[i] = (double *) malloc(ny * sizeof(double));
      if (im[i] == NULL) {
	fprintf (stderr, "ERROR en toImage, im[%d] = NULL.\n", i);
	exit (1);
      }
      for (j = 0; j < ny; j++) {
	im[i][j] = NULL_PIX;
	mask[i][j] = -1;
      }
    }
  }
  
  for (n = malla->aristas->primero; n != NULL; n = n->sgte) {
    a = (AristaVoronoi *) n->info;
    if ((fabs(a->ptoIni->x - a->ptoFin->x) < 0.5/nx) &&
	(fabs(a->ptoIni->y - a->ptoFin->y) < 0.5/ny)) {
      continue;
    }
    m = (a->ptoIni->y - a->ptoFin->y) * ny / 
      ((a->ptoIni->x - a->ptoFin->x) * nx);
    if (a->ptoIni->y < a->ptoFin->y) {
      polD = a->poliDer;
      polI = a->poliIzq;
    }
    else if (a->ptoIni->y > a->ptoFin->y) {
      polI = a->poliDer;
      polD = a->poliIzq;
    }
    else { // ==
      if (a->ptoIni->x < a->ptoFin->x) {
	polI = a->poliDer;
	polD = a->poliIzq;	
      }
      else {
	polD = a->poliDer;
	polI = a->poliIzq;	
      }
    }
    if (m > 1) {
      linea3(im, nx, ny, a->ptoIni, a->ptoFin, polD, polI, mask);
    }
    else if (m > 0) {
      linea1(im, nx, ny, a->ptoFin, a->ptoIni, polD, polI, mask);
    }
    else if (m > -1) {
      linea2(im, nx, ny, a->ptoIni, a->ptoFin, polD, polI, mask);
    }
    else {
      linea4(im, nx, ny, a->ptoIni, a->ptoFin, polD, polI, mask);
    }
  }

  raster (im, mask, nx, ny, malla);

  return im;  
}

double **toImageSinRaster(MallaVoronoi *malla, int nx, int ny, int **mask) {
  double **im = (double **) malloc(nx * sizeof(double *));
  int i, j;
  double m;
  NodoLista *n;
  AristaVoronoi *a;
  PoligonoVoronoi *polD, *polI;

  if (im == NULL) {
    fprintf (stderr, "ERROR en toImage, im = NULL.\n");
    exit (1);
  }

  if (mask == NULL) {
    for(i = 0; i < nx; i++) {
      im[i] = (double *) malloc(ny * sizeof(double));
      if (im[i] == NULL) {
	fprintf (stderr, "ERROR en toImage, im[%d] = NULL.\n", i);
	exit (1);
      }
      for (j = 0; j < ny; j++) {
	im[i][j] = NULL_PIX;
      }
    }
  }
  else {
    for(i = 0; i < nx; i++) {
      im[i] = (double *) malloc(ny * sizeof(double));
      if (im[i] == NULL) {
	fprintf (stderr, "ERROR en toImage, im[%d] = NULL.\n", i);
	exit (1);
      }      
      for (j = 0; j < ny; j++) {
	im[i][j] = NULL_PIX;
	mask[i][j] = -1;
      }
    }
  }
  
  for (n = malla->aristas->primero; n != NULL; n = n->sgte) {
    a = (AristaVoronoi *) n->info;
    if ((fabs(a->ptoIni->x - a->ptoFin->x) < 0.5/nx) &&
	(fabs(a->ptoIni->y - a->ptoFin->y) < 0.5/ny)) {
      continue;
    }
    m = (a->ptoIni->y - a->ptoFin->y) * ny / 
      ((a->ptoIni->x - a->ptoFin->x) * nx);
    if (a->ptoIni->y < a->ptoFin->y) {
      polD = a->poliDer;
      polI = a->poliIzq;
    }
    else if (a->ptoIni->y > a->ptoFin->y) {
      polI = a->poliDer;
      polD = a->poliIzq;
    }
    else { // ==
      if (a->ptoIni->x < a->ptoFin->x) {
	polI = a->poliDer;
	polD = a->poliIzq;	
      }
      else {
	polD = a->poliDer;
	polI = a->poliIzq;	
      }
    }
    if (m > 1) {
      linea3(im, nx, ny, a->ptoIni, a->ptoFin, polD, polI, mask);
    }
    else if (m > 0) {
      linea1(im, nx, ny, a->ptoFin, a->ptoIni, polD, polI, mask);
    }
    else if (m > -1) {
      linea2(im, nx, ny, a->ptoIni, a->ptoFin, polD, polI, mask);
    }
    else {
      linea4(im, nx, ny, a->ptoIni, a->ptoFin, polD, polI, mask);
    }
  }

  return im;  
}

void linea1(double **im, int nx, int ny,
	    PuntoVoronoi *p1, 
	    PuntoVoronoi *p2,
	    PoligonoVoronoi *polD,
	    PoligonoVoronoi *polI,
	    int **mask) {
  double dx, dy, d, incrE, incrNE, ay, x, y;
  int ifin, i, j;
  PuntoVoronoi *pini, *pfin, *qini, *qfin;

  if(p1->x < p2->x) {
    pini = p1;
    pfin = p2;
  }
  else {
    pini = p2;
    pfin = p1;
  }

  dy = (pfin->y - pini->y) * (ny - 1);
  dx = (pfin->x - pini->x) * (nx - 1);
  incrE = dy * 2.0;
  incrNE = (dy - dx) * 2.0;

  qini = (PuntoVoronoi *) malloc (sizeof (PuntoVoronoi));
  if (qini == NULL) {
    fprintf (stderr, "ERROR en linea1, qini = NULL.\n");
    exit (1);
  }

  qfin = (PuntoVoronoi *) malloc (sizeof (PuntoVoronoi));
  if (qfin == NULL) {
    fprintf (stderr, "ERROR en linea1, qfin = NULL.\n");
    exit (1);
  }

  if (!interseccionCuadrado (pini, pfin, qini, qfin, 1)) {
    free(qini);
    free(qfin);
    return;
  }

  ifin = (int) round(qfin->x * (nx - 1) + 1);
  x = round(qini->x * (nx - 1) - 1);
  y = (x - qini->x * (nx - 1)) * dy / dx + qini->y * (ny - 1) + 0.5;
  i = (int) x;
  j = (int) round(y);
  ay = (j + 0.5) - y;
  d = 2.0 * (dy - dx * ay);

  for (; i <= ifin + 1 && i <= nx && j <= ny; i++) {
    /*
    y = (x - qini->x * (nx - 1)) * dy / dx + qini->y * (ny - 1) + 0.5;
    if (i != round(x) || j != round(y)) {
      fprintf (stderr, "ERROR en linea1: (x, y) = (%g, %g)\n", x, y);
      fprintf (stderr, "                 (i, j) = (%d, %d)\n", i, j);
      fprintf (stderr, "                 p1, p2 = (%g, %g), (%g, %g)\n",
	       p1->x, p1->y, p2->x, p2->y);
      f = 2.0 * (dy * (i + 1.0) - dx * (j + 0.5) + 
		 (qini->y * (ny - 1) + 0.5) * dx - (qini->x * (nx - 1))* dy);
      fprintf (stderr, "              d      = %g\n", d);
      fprintf (stderr, "              f      = %g\n", f);
      exit (1);
    }
    x += 1.0;
    */
    if (i < nx && i >= 0) {
      if(j < ny && j >= 0) {
	if (dentro (polI, i / (nx - 1.0), j / (ny - 1.0))) {
	  im[i][j] = polI->valor;
	  if (mask) {
	    mask[i][j] = polI->id;
	  }
	}
      }
      if(j - 1 < ny && j - 1 >= 0) {
	if (dentro (polD, i / (nx - 1.0), (j - 1.0) / (ny - 1.0))) {
	  im[i][j - 1] = polD->valor;
	  if (mask) {
	    mask[i][j - 1] = polD->id;
	  }
	}
      }
    }
    if(d < 0) {
      d += incrE;
    }
    else {
      d += incrNE;
      j++;
    }
  }
  free (qini);
  free (qfin);
}

void linea2(double **im, int nx, int ny,
	     PuntoVoronoi *p1, 
	     PuntoVoronoi *p2,
	     PoligonoVoronoi *polD,
	     PoligonoVoronoi *polI,
	    int **mask) {
  double dx, dy, d, incrE, incrNE, ay, x, y;
  int ifin, i, j;
  PuntoVoronoi *pini, *pfin, *qini, *qfin;

  if(p1->x < p2->x) {
    pini = p1;
    pfin = p2;
  }
  else {
    pini = p2;
    pfin = p1;
  }

  dy = (pfin->y - pini->y) * (ny - 1);
  dx = (pfin->x - pini->x) * (nx - 1);
  incrE = dy * 2.0;
  incrNE = (dy + dx) * 2.0;

  qini = (PuntoVoronoi *) malloc (sizeof (PuntoVoronoi));
  if (qini == NULL) {
    fprintf (stderr, "ERROR en linea1, qini = NULL.\n");
    exit (1);
  }
  qfin = (PuntoVoronoi *) malloc (sizeof (PuntoVoronoi));
  if (qfin == NULL) {
    fprintf (stderr, "ERROR en linea1, qfin = NULL.\n");
    exit (1);
  }
  if (!interseccionCuadrado (pini, pfin, qini, qfin, 2)) {
    free(qini);
    free(qfin);
    return;
  }

  ifin = (int) round(qfin->x * (nx - 1) + 1);

  x = round(qini->x * (nx - 1) - 1);
  y = (x - qini->x * (nx - 1)) * dy / dx + qini->y * (ny - 1) + 0.5;
  i = (int) x;
  j = (int) round(y);
  ay = y + 0.5 - j;
  d = 2.0 * (dy  + dx * ay);

  for (; i <= ifin && i <= nx && j >= -1; i++) {
    /*    
    y = (x - qini->x * (nx - 1)) * dy / dx + qini->y * (ny - 1) + 0.5;
    if (i != round(x) || j != round(y)) {
      fprintf (stderr, "ERROR en linea2: (x, y) = (%g, %g)\n", x, y);
      fprintf (stderr, "                 (i, j) = (%d, %d)\n", i, j);
      fprintf (stderr, "                 p1, p2 = (%g, %g), (%g, %g)\n",
	       p1->x, p1->y, p2->x, p2->y);
      exit (1);
    }
    x += 1.0;
    */
    if (i < nx && i >= 0) {
      if(j < ny && j >= 0) {
	if (dentro (polD, i / (nx - 1.0), j / (ny - 1.0))) {
	  im[i][j] = polD->valor;
	  if (mask) {
	    mask[i][j] = polD->id;
	  }
	}
      }
      if(j - 1 < ny && j - 1 >= 0) {
	if (dentro (polI, i / (nx - 1.0), (j - 1.0) / (ny - 1.0))) {
	  im[i][j - 1] = polI->valor;
	  if (mask) {
	    mask[i][j - 1] = polI->id;
	  }
	}
      }
    }

    if(d < 0) {
      d += incrNE;
      j--;
    }
    else {
      d += incrE;
    }
  }
  free (qini);
  free (qfin);
}

void linea3(double **im, int nx, int ny,
	    PuntoVoronoi *p1, 
	    PuntoVoronoi *p2,
	    PoligonoVoronoi *polD,
	    PoligonoVoronoi *polI,
	    int **mask) {
  double dx, dy, d, incrE, incrNE, ax, x, y;
  int jfin, i, j;
  PuntoVoronoi *pini, *pfin, *qini, *qfin;

  if(p1->y < p2->y) {
    pini = p1;
    pfin = p2;
  }
  else {
    pini = p2;
    pfin = p1;
  }

  dy = (pfin->y - pini->y) * (ny - 1);
  dx = (pfin->x - pini->x) * (nx - 1);
  incrE = - dx * 2.0;
  incrNE = (dy - dx) * 2.0;

  qini = (PuntoVoronoi *) malloc (sizeof (PuntoVoronoi));
  if (qini == NULL) {
    fprintf (stderr, "ERROR en linea1, qini = NULL.\n");
    exit (1);
  }
  qfin = (PuntoVoronoi *) malloc (sizeof (PuntoVoronoi));
  if (qfin == NULL) {
    fprintf (stderr, "ERROR en linea1, qfin = NULL.\n");
    exit (1);
  }
  if (!interseccionCuadrado (pini, pfin, qini, qfin, 3)) {
    free(qini);
    free(qfin);
    return;
  }

  jfin = (int) round (qfin->y * (ny - 1));

  y = round (qini->y * (ny - 1) - 1);
  x = (y - qini->y * (ny - 1)) * dx / dy + qini->x * (nx - 1) + 0.5;
  i = (int) round (x);
  j = (int) y;
  ax = i + 0.5 - x;
  d = 2.0 * (dy * ax - dx);

  for (; j <= jfin + 1 && j <= ny && i <= nx; j++) {
    /*      
    x = (y - qini->y * (ny - 1)) * dx / dy + qini->x * (nx - 1) + 0.5;
    if (i != round(x) || j != round(y)) {
      fprintf (stderr, "ERROR en linea3: (x, y) = (%g, %g)\n", x, y);
      fprintf (stderr, "                 (i, j) = (%d, %d)\n", i, j);
      fprintf (stderr, "                 p1, p2 = (%g, %g), (%g, %g)\n",
	       p1->x, p1->y, p2->x, p2->y);
      exit (1);
    }
    y += 1.0;
    */
    if(j < ny && j >= 0) {
      if (i < nx && i >= 0) {
	if (dentro (polD, i / (nx - 1.0), j / (ny - 1.0))) {
	  im[i][j] = polD->valor;
	  if (mask) {
	    mask[i][j] = polD->id;
	  }
	}
      }
      if(i - 1 < nx && i - 1 >= 0) {
	if (dentro (polI, (i - 1.0) / (nx - 1.0), j / (ny - 1.0))) {
	  im[i - 1][j] = polI->valor;
	  if (mask) {
	    mask[i - 1][j] = polI->id;
	  }
	}
      }
    }
    if(d < 0) {
      d += incrNE;
      i++;
    }
    else {
      d += incrE;
    }
  }
  if (i == nx && j < ny && j >= 0) {
    im[i - 1][j] = polI->valor;
    if (mask) {
      mask[i - 1][j] = polI->id;
    }
  }
  free (qini);
  free (qfin);
}

void linea4(double **im, int nx, int ny,
	    PuntoVoronoi *p1, 
	    PuntoVoronoi *p2,
	    PoligonoVoronoi *polD,
	    PoligonoVoronoi *polI,
	    int **mask) {
  double dx, dy, d, incrE, incrNE, ax, x, y;
  int jfin, i, j;
  PuntoVoronoi *pini, *pfin, *qini, *qfin;

  if(p1->y > p2->y) {
    pini = p1;
    pfin = p2;
  }
  else {
    pini = p2;
    pfin = p1;
  }
  
  dy = (pfin->y - pini->y) * (ny - 1);
  dx = (pfin->x - pini->x) * (nx - 1);
  incrE = dx * 2.0;
  incrNE = (dy + dx) * 2.0;

  qini = (PuntoVoronoi *) malloc (sizeof (PuntoVoronoi));
  if (qini == NULL) {
    fprintf (stderr, "ERROR en linea1, qini = NULL.\n");
    exit (1);
  }
  qfin = (PuntoVoronoi *) malloc (sizeof (PuntoVoronoi));
  if (qfin == NULL) {
    fprintf (stderr, "ERROR en linea1, qfin = NULL.\n");
    exit (1);
  }
  if (!interseccionCuadrado (pini, pfin, qini, qfin, 4)) {
    free(qini);
    free(qfin);
    return;
  }
  
  jfin = (int) round (qfin->y * (ny - 1));
  
  y = round (qini->y * (ny - 1) + 1);
  x = (y - qini->y * (ny - 1)) * dx / dy + qini->x * (nx - 1) + 0.5;
  i = (int) round(x);
  j = (int) y;
  ax = i + 0.5 - x;
  //ax = pini->x * nx - (i - 0.5);
  //ay = j + 1.0 - pini->y * ny;  
  d = 2.0 * (dx + dy * ax);
  //d = - 2.0 * (dx * (ay + 0.5) + dy * ax);

  for (; j >= jfin - 1 && j >= 0 && i <= nx; j--) {
    /*
    x = (y - qini->y * (ny - 1)) * dx / dy + qini->x * (nx - 1) + 0.5;
    if (i != round (x) || j != round (y)) {
      fprintf (stderr, "ERROR en linea4: (x, y) = (%g, %g)\n", x, y);
      fprintf (stderr, "                 (i, j) = (%d, %d)\n", i, j);
      fprintf (stderr, "                 p1, p2 = (%g, %g), (%g, %g)\n",
	       p1->x, p1->y, p2->x, p2->y);
      f = 2.0 * (dy * (i + 1.0) - dx * (j - 0.5) + 
		 pini->y * dx * (ny - 1) - (pini->x * (nx - 1) + 0.5)* dy);
      fprintf (stderr, "              d      = %g\n", d);
      fprintf (stderr, "              f      = %g\n", f);
      
      exit (1);
    }
    
    y -= 1.0;
    */
    if(j < ny && j >= 0) {
      if (i < nx && i >= 0) {
	if (dentro (polD, i / (nx - 1.0), j / (ny - 1.0))) {
	  im[i][j] = polD->valor;
	  if (mask) {
	    mask[i][j] = polD->id;
	  }
	}
      }
      if(i - 1 < nx && i - 1 >= 0) {
	if (dentro (polI, (i - 1.0) / (nx - 1.0), j / (ny - 1.0))) {
	  im[i - 1][j] = polI->valor;
	  if (mask) {
	    mask[i - 1][j] = polI->id;
	  }
	}
      }
    }
    if(d < 0) {
      d += incrE;
    }
    else {
      d += incrNE;
      i++;
    }
  }

  free (qini);
  free (qfin);
}

/* Calcula las intersecciones del cuadrado (0,0) (0.1) con la recta formada 
 * por pini y pfin dejando los ptos finales e iniciales en qfin y qini
 * respectivamente. tipoLinea es el tipo de linea definido por su pendiente.
 * En caso de no haber interseccion retorna FALSE.
 */
int interseccionCuadrado (PuntoVoronoi *pini, PuntoVoronoi *pfin,
			  PuntoVoronoi *qini, PuntoVoronoi *qfin,
			  int tipoLinea) {
  double m;

  if (pini->x > 1 || pfin->x < 0) {
    return FALSE;
  }
  if ((pini->y < 0 && pfin->y < 0) || (pini->y > 1 && pfin->y > 1)) {
    return FALSE;
  }
  if (pini->y == pfin->y) { // Linea Horizontal
    if (pini->x < 0) {
      qini->x = 0;
    }
    else {
      qini->x = pini->x;    
    }
    qini->y = pini->y;

    if (pfin->x > 1) {
      qfin->x = 1;
    }
    else {
      qfin->x = pfin->x;
    }
    qfin->y = pini->y;
    return TRUE;
  }
  
  if (tipoLinea  == 1 || tipoLinea == 3) { // Pendiente > 0
    if (pini->x == pfin->x) { // Linea vertical
      qini->x = pini->x;
      if (pini->y < 0) {
	qini->y = 0;
      }
      else {
	qini->y = pini->y;
      }
      qfin->x = pfin->x;
      if (pfin->y > 1) {
	qfin->y = 1;
      }
      else {
	qfin->y = pfin->y;
      }
      return TRUE;
    }

    m = (pfin->y - pini->y) / (pfin->x - pini->x);
    
    qini->x = (pini->x < 0)? 0 : pini->x;
    qfin->x = (pfin->x > 1)? 1 : pfin->x;
    
    qini->y = m * (qini->x - pini->x) + pini->y;
    qfin->y = m * (qfin->x - pini->x) + pini->y;
    
    if (qini->y < 0) { 
      qini->y = 0;
      qini->x = (qini->y - pini->y) / m + pini->x;
    }
    else if (qini->y > 1) {
      return FALSE;
      //qini->y = 1;
      //pini->x = dx/dy * (qini->y - pini->y) + pini->x;
    }
    if (qini->x > 1) {
      return FALSE;
    }
    
    if (qfin->y < 0) {
      return FALSE;
      //qfin->y = 0;
      //qfin->x = dx/dy * (qfin->y - pini->y) + pini->x;
    }
    else if (qfin->y > 1) {
      qfin->y = 1;
      qfin->x = (qfin->y - pini->y) / m + pini->x;
    }
  }
  else if (tipoLinea  == 2 || tipoLinea == 4) { // Pendiente < 0
    if (pini->x == pfin->x) { // Linea vertical
      qini->x = pini->x;
      if (pini->y >= 1) {
	qini->y = 1;
      }
      else {
	qini->y = pini->y;
      }
      qfin->x = pfin->x;
      if (pfin->y < 0) {
	qfin->y = 0;
      }
      else {
	qfin->y = pfin->y;
      }
      return TRUE;
    }

    m = (pfin->y - pini->y) / (pfin->x - pini->x);
    
    qini->x = (pini->x < 0)? 0 : pini->x;
    qfin->x = (pfin->x > 1)? 1 : pfin->x;
    
    qini->y = m * (qini->x - pini->x) + pini->y;
    qfin->y = m * (qfin->x - pini->x) + pini->y;
    
    if (qini->y < 0) { 
      return FALSE;
      //qini->y = 0;
      //qini->x = (qini->y - pini->y) / m + pini->x;
    }
    else if (qini->y > 1) {
      //return FALSE;
      qini->y = 1;
      qini->x = (qini->y - pini->y) / m + pini->x;
    }
    if (qini->x > 1) {
      return FALSE;
    }
    
    if (qfin->y < 0) {
      //return FALSE;
      qfin->y = 0;
      qfin->x = (qfin->y - pini->y) / m + pini->x;
    }
    else if (qfin->y > 1) {
      return FALSE;
      //qfin->y = 1;
      //qfin->x = (qfin->y - pini->y) / m + pini->x;
    }
  }
  else {
    return FALSE;
  }
  return TRUE;
}

void raster (double **im, int **mask, int nx, int ny, MallaVoronoi *m) {
  int i, j;
  double valorIm;
  int valorMask = -1;

  if(bordeV(im, mask, 0, ny)) {
    //Comenzamos raster horizontal

    for(j = 0; j < ny; j++) {
      valorIm = im[0][j];
      if (mask) {
	valorMask = mask[0][j];
      }
      for(i = 0; i < nx; i++) {
	if (im[i][j] != NULL_PIX) {
	  valorIm = im[i][j];
	  if (mask) {
	    valorMask = mask[i][j];
	  }
	}
	else {
	  im[i][j] = valorIm;
	  if (mask) {
	    /*
	    double x, y;
	    PoligonoVoronoi *pol;
	    x = (double) i / (nx - 1);
	    y = (double) j / (ny - 1);
	    pol = encontrarPoligono (m, x, y);
	    if (pol->id != valorMask) {
	      fprintf (stderr, "ERROR en raster (%d, %d)\n", i, j);
	      fprintf (stderr, "      pol->id = %d\n", pol->id);
	      fprintf (stderr, "      valorMask = %d\n", valorMask);
	      fprintf (stderr, "      NULL_PIX = %g\n", NULL_PIX);
	      exit (1);
	      }*/
	    mask[i][j] = valorMask;
	  }
	}
      }
    }
  }

  else if (bordeV(im, mask, nx - 1, ny)) {
    //Comenzamos raster horizontal inverso
    //printf ("Borde V derecha\n");

    for(j = 0; j < ny; j++) {
      valorIm = im[0][j];
      if (mask) {
	valorMask = mask[0][j];
      }
      for(i = nx - 1; i >= 0; i--) {
	if (im[i][j] != NULL_PIX) {
	  valorIm = im[i][j];
	  if (mask) {
	    valorMask = mask[i][j];
	  }
	}
	else {
	  im[i][j] = valorIm;
	  if (mask) {
	    mask[i][j] = valorMask;
	  }
	}
      }
    }
  }

  else if (bordeH(im, mask, nx, 0)) {
    //Comenzamos raster vertical
    //printf ("Borde H arriba\n");
    
    for(i = 0; i < nx; i++) {
      valorIm = im[i][0];
      if (mask) {
	valorMask = mask[i][0];
      }
      for(j = 0; j < ny; j++) {
	if (im[i][j] != NULL_PIX) {
	  valorIm = im[i][j];
	  if (mask) {
	    valorMask = mask[i][j];
	  }
	}
	else {
	  im[i][j] = valorIm;
	  if (mask) {
	    mask[i][j] = valorMask;
	  }
	}
      }
    }
  }  
}

int bordeV (double **im, int **mask, int i, int ny) {
  int j1, j2;
  double valorIm;
  int valorMask = -1;
  
  for (j1 = 0; j1 < ny && im[i][j1] == NULL_PIX; j1++);

  if (j1 < ny) {
    valorIm = im[i][j1];
    if (mask) {
      valorMask = mask[i][j1];
    }
    
    for (j2 = j1 - 1; j2 >= 0; j2--){
      im[i][j2] = valorIm;
      if (mask) {
	mask[i][j2] = valorMask;
      }
    }
    for(++j1; j1 < ny; j1++) {
      if (im[i][j1] != NULL_PIX) {
	valorIm = im[i][j1];
	if (mask) {
	  valorMask = mask[i][j1];
	}
      }
      else {
	im[i][j1] = valorIm;
	if (mask) {
	  mask[i][j1] = valorMask;
	}
      }
    }
    return TRUE;
  }
  return FALSE;
}

int bordeH (double **im, int **mask, int nx, int j) {
  int i1, i2;
  double valorIm;
  int valorMask = -1;
  
  for (i1 = 0; i1 < nx && im[i1][j] == NULL_PIX; i1++);
  
  if (i1 < nx) {
    valorIm = im[i1][j];
    if (mask) {
      valorMask = mask[i1][j];
    }
    for (i2 = i1 - 1; i2 >= 0; i2--){
      im[i2][j] = valorIm;
      if (mask) {
	mask[i2][j] = valorMask;
      }
    }
    for(++i1; i1 < nx; i1++) {
      if (im[i1][j] != NULL_PIX) {
	valorIm = im[i1][j];
	if (mask) {
	  valorMask = mask[i1][j];
	}
      }
      else {
	im[i1][j] = valorIm;
	if (mask) {
	  mask[i1][j] = valorMask;
	}
      }
    }
    return TRUE;
  }
  return FALSE;
}

int comprobarImagen (MallaVoronoi *m, int **mask, int nx, int ny) {
  int i, j, ret = TRUE;
  double x, y;
  PoligonoVoronoi * pol;

  for (i = 0; i < nx; i++) {
    x = (double) i / (nx - 1);
    for (j = 0; j < ny; j++) {
      y = (double) j / (ny - 1);
      pol = encontrarPoligono (m, x, y);
      if (pol->id != mask [i][j]){
	ret = FALSE;
	fprintf (stderr, "ERROR al comprobarImagen\n");
	fprintf (stderr, "      en pixel (%d, %d)\n", i, j);
	imprimirMallaArchivo (m, "ERRORcomprobarImagen.dat");
      }
    }
  }
  return ret;
}
/*
void submkim(float *xcel, float *ycel, float *icel, int ncel, int nx, int ny,
	     float im[nx][ny]) 
{
  MallaVoronoi *mesh = newMallaVoronoi();
  int i, **mask;
  double **im2;
  int j;
   
  for (i = 0; i < ncel; i++)
    {
      insertarSitio (mesh, xcel[i], ycel[i], icel[i]);
    }


  mask = (int **) malloc(nx * sizeof(int *));
  for (i = 0; i < nx; i++) {
    mask[i] = (int *) malloc (ny * sizeof(int));
  }

  im2 = toImage (mesh, nx, ny, mask);

  imprimirAristasMallaArchivo(mesh, "aristas.dat", nx, ny);  

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      im[i][j] = im2[i][j];
    }
    free (im2[i]);
  }
  free (im2);
}
*/
