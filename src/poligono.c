#include "poligono.h"
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

#define DEBUG 0
#define DEBUG2 0 // imprime las aristas que se van sumando

PoligonoVoronoi *newPoligonoVoronoi(double x, double y, double valor, int id) {
  PoligonoVoronoi *p = 
    (PoligonoVoronoi *) malloc(sizeof(PoligonoVoronoi));
  p->x = x;
  p->y = y;
  p->valor = valor;
  p->id = id;
  return p;
}

int dentro (PoligonoVoronoi *pol, double x, double y) {
  AristaVoronoi *a, *aIni;
  int der, izq;
  double x1;

  der = FALSE;
  izq = FALSE;
  aIni = pol->a;
  a = aIni;
  do {
    if ((a->ptoIni->y <= y && y <= a->ptoFin->y) ||
	(a->ptoFin->y <= y && y <= a->ptoIni->y)) {
      // Arista a la altura del punto.
      
      if (a->ptoIni->y == a->ptoFin->y) {
	return TRUE;
      }
      x1 = ((y - a->ptoIni->y) * (a->ptoIni->x - a->ptoFin->x) / 
	    (a->ptoIni->y - a->ptoFin->y) + a->ptoIni->x);

      if (x <= x1) {
	if (izq) {
	  return TRUE;
	  }
	der = TRUE;
      }
      else if (x1 <= x) {
	if (der) {
	  return TRUE;
	}
	izq = TRUE;
      }
    }
    if (a->poliDer == pol) {
      a = a->ccwSucc;
    }
    else {
      a = a->ccwPred;  
    }
  }
  while (a != aIni);
  return FALSE;
}

double area (PoligonoVoronoi *pol) {
  double suma = 0;
  AristaVoronoi *a;

  a = pol->a;
  do {
    if (a->poliDer == pol) {
      suma += a->ptoFin->x * a->ptoIni->y - a->ptoIni->x * a->ptoFin->y;
      a = a->cwPred;
    }
    else if (a->poliIzq == pol) {
      suma += a->ptoIni->x * a->ptoFin->y - a->ptoFin->x * a->ptoIni->y;
      a = a->cwSucc;
    }
    else {
      printf ("ERROR en area: Mal recorrido de aristas.\n");
      exit (1);
    }
  }
  while (a != pol->a);

  return suma / 2;
}
  
double areaCuadrado (PoligonoVoronoi *pol){
  double suma = 0;
  AristaVoronoi *a, *a_old;
  PuntoVoronoi *p1, *p2, *q1, *q2;

  a = pol->a;

  // para comenzar dentro del cuadrado...
  if (DEBUG) {
    printf("dentro?\n");
  }
  do {
    if (a->poliDer == pol) {
      a = a->cwPred;
    }
    else if (a->poliIzq == pol) {
      a = a->cwSucc;
    }
    if (a->poliDer == pol) {
      p1 = a->ptoFin;
      p2 = a->ptoIni;
    }
    else if (a->poliIzq == pol) {
      p1 = a->ptoIni;
      p2 = a->ptoFin;
    }
    else {
      printf ("ERROR en areaCuadrado, arista no corresponde a pol.\n");
      exit (1);
    }
    if (DEBUG) {
      printf ("%d ", dentroCuadrado(a));
      imprimirArista(a);
      }
  } while (!dentroCuadrado(a) && a != pol->a);

  if (!dentroCuadrado(a)) {
    return FALSE;
  }

  a_old = a;
  q1 = newPuntoVoronoi(p1->x, p1->y);
  q2 = newPuntoVoronoi(p2->x, p2->y);

  if (DEBUG) {
    printf("dentro\n");
  }

  do {
    if (bordeCuadrado (p1, p2, q1, q2) && 
	(p2->x > 1 || p2->x < 0 || p2->y > 1 || p2->y < 0)) {

      if (p1->x > p2->x) { // bordeCuadrado los dio vuelta
	suma += q1->y * q2->x - q2->y * q1->x;
	if (DEBUG2) printf("sumando0.1 (%g, %g),(%g, %g)\n", q2->x, q2->y, q1->x, q1->y);
      }
      else {
	suma += q1->x * q2->y - q2->x * q1->y;
	if (DEBUG2) printf("sumando0.2 (%g, %g),(%g, %g)\n", q1->x, q1->y, q2->x, q2->y);
      }

      suma += sumarBordes(&a, p1, p2, q1, q2, pol);
    }
    else {
      if (p1->x > p2->x) { // bordeCuadrado los dio vuelta
	suma += q1->y * q2->x - q2->y * q1->x;
	if (DEBUG2) printf("sumando3 (%g, %g),(%g, %g)\n", q2->x, q2->y, q1->x, q1->y);
      }
      else {
	suma += q1->x * q2->y - q2->x * q1->y;
	if (DEBUG2) printf("sumando4 (%g, %g),(%g, %g)\n", q1->x, q1->y, q2->x, q2->y);
      }
      if (a->poliDer == pol) {
	a = a->cwPred;
      }
      else if (a->poliIzq == pol) {
	a = a->cwSucc;
      }
    }
    if (a->poliDer == pol) {
      p1 = a->ptoFin;
      p2 = a->ptoIni;
    }
    else if (a->poliIzq == pol) {
      p1 = a->ptoIni;
      p2 = a->ptoFin;
    }
    else {
      printf ("ERROR en area: Mal recorrido de aristas.\n");
      exit (1);
    }
  } while (a != a_old);

  free (q1);
  free (q2);
  return suma / 2;
}

/* Calcula las intersecciones del cuadrado (0,0) (0.1) con la recta
 * formada por pini y pfin dejando los ptos finales e iniciales en
 * qfin y qini respectivamente. tipoLinea es el tipo de linea definido
 * por su pendiente.  En caso de no haber interseccion retorna FALSE.
 */
int bordeCuadrado (PuntoVoronoi *pini1, PuntoVoronoi *pfin1,
		   PuntoVoronoi *qini, PuntoVoronoi *qfin) {
  double m;
  int tipoLinea, cambiado = FALSE;
  PuntoVoronoi *pini, *pfin;

  m = (pfin1->y - pini1->y) / (pfin1->x - pini1->x);

  if (m > 1) {
    tipoLinea = 3;
  }
  else if (m > 0) {
    tipoLinea = 1;
  }
  else if (m > -1) {
    tipoLinea = 2;
  }
  else {
    tipoLinea = 4;
  }

  if (DEBUG) printf ("tipo Linea = %d\n", tipoLinea);

  // ordenamos segun y
  if (pini1->y > pfin1->y) {
    pini = pfin1;
    pfin = pini1;
  }
  else {
    pini = pini1;
    pfin = pfin1;
  }
  if (pini->y > 1 || pfin->y < 0) {
    if (pini1->x > pfin1->x) {
      pini = pfin1;
      pfin = pini1;
    }
    else {
      pini = pini1;
    pfin = pfin1;
    }
    
    qini->x = pini->x;
    qini->y = pini->y;
    qfin->x = pfin->x;
    qfin->y = pfin->y;

    return FALSE;
  }
  
  //ordenamos segun y
  if (pini1->x > pfin1->x) {
    pini = pfin1;
    pfin = pini1;
  }
  else {
    pini = pini1;
    pfin = pfin1;
  }

  qini->x = pini->x;
  qini->y = pini->y;
  qfin->x = pfin->x;
  qfin->y = pfin->y;

  if (pini->x > 1 || pfin->x < 0) {
    return FALSE;
  }

  if (pini->y == pfin->y) { // Linea Horizontal
    if (pini->x < 0) {
      qini->x = 0;
      cambiado = TRUE;
    }
    if (qfin->x > 1) {
      qfin->x = 1;
      cambiado = TRUE;
    }
    return cambiado;
  }
  
  if (tipoLinea  == 1 || tipoLinea == 3) { // Pendiente > 0
    if (pini->x == pfin->x) { // Linea vertical
      if (DEBUG) printf ("vertical\n");
      if (pini->y < 0) {
	qini->y = 0;
	cambiado = TRUE;
      }
      if (pfin->y > 1) {
	qfin->y = 1;
	cambiado = TRUE;
      }
      return cambiado;
    }

    if (pini->x < 0) {
      qini->x = 0;
      cambiado = TRUE;
    }    
    if (pfin->x > 1) {
      qfin->x = 1;
      cambiado = TRUE;
    }    
    
    qini->y = m * (qini->x - pini->x) + pini->y;
    qfin->y = m * (qfin->x - pini->x) + pini->y;
    
    if (qini->y < 0) { 
      qini->y = 0;
      qini->x = (qini->y - pini->y) / m + pini->x;
      cambiado = TRUE;
    }
    else if (qini->y > 1) {
      return FALSE;
    }
    if (qini->x > 1) {
      return FALSE;
    }
    else if (qfin->y < 0) {
      return FALSE;
    }
    else {
      if (qfin->y > 1) {
	qfin->y = 1;
	qfin->x = (qfin->y - pini->y) / m + pini->x;
	cambiado = TRUE;
      }
      return cambiado;
    }
  }
  else if (tipoLinea  == 2 || tipoLinea == 4) { // Pendiente < 0
    if (pini->x == pfin->x) { // Linea vertical
      if (DEBUG) printf ("vertical\n");
      if (pini->y >= 1) {
	qini->y = 1;
	cambiado = TRUE;
      }
      if (pfin->y < 0) {
	qfin->y = 0;
	cambiado = TRUE;
      }
      return cambiado;
    }
    
    if (pini->x < 0) {
      qini->x = 0;
      cambiado = TRUE;
    }
    if (pfin->x > 1) {
      qfin->x = 1;
      cambiado = TRUE;
    }
    
    qini->y = m * (qini->x - pini->x) + pini->y;
    qfin->y = m * (qfin->x - pini->x) + pini->y;
    
    if (qini->y < 0) { 
      return FALSE;
    }
    else if (qini->y > 1) {
      qini->y = 1;
      qini->x = (qini->y - pini->y) / m + pini->x;
      cambiado = TRUE;
    }
    if (qini->x > 1) {
      return FALSE;
    }
    
    if (qfin->y < 0) {
      qfin->y = 0;
      qfin->x = (qfin->y - pini->y) / m + pini->x;
      return TRUE;
    }
    else if (qfin->y > 1) {
      return FALSE;
    }
  }
  else {
    return cambiado;
  }
  return cambiado;
}

int dentroCuadrado(AristaVoronoi *a) {
  double dx, dy, f1, f2, x1, y1, x2, y2;
  PuntoVoronoi *p1, *p2;

  //ordenamos segun x
  if (a->ptoIni->x < a->ptoFin->x) {
    p1 = a->ptoIni;
    p2 = a->ptoFin;
  }
  else {
    p1 = a->ptoFin;    
    p2 = a->ptoIni;
  }

  if (p1->x > 1 || p2->x < 0) {
    return FALSE;
  }
  
  //ordenamos segun y
  if (a->ptoIni->y < a->ptoFin->y) {
    p1 = a->ptoIni;
    p2 = a->ptoFin;
  }
  else {
    p1 = a->ptoFin;    
    p2 = a->ptoIni;
  }

  if (p1->y > 1 || p2->y < 0) {
    return FALSE;
  }

  dx = p1->x - p2->x;
  dy = p1->y - p2->y;
  if (dy / dx > 0) {
    x1 = 1;
    y1 = 0;
    x2 = 0;
    y2 = 1;
  }
  else {
    x1 = 0;
    y1 = 0;
    x2 = 1;
    y2 = 1;
  }
  f1 = dy * (x1 - a->ptoIni->x) - dx * (y1 - a->ptoIni->y);
  f2 = dy * (x2 - a->ptoIni->x) - dx * (y2 - a->ptoIni->y);

  return (f1 * f2 <= 0);
}

int tope (PuntoVoronoi *p) {
  if (p->x == 0) {
    return BORDE_IZQ;
  } 
  else if (p->x == 1) {
    return BORDE_DER;
  }
  else if (p->y == 0) {
    return BORDE_INF;
  }
  else if (p->y == 1) {
    return BORDE_SUP;
  }
  return 0;
}

/* Suma los bordes de comenzando por 'a' hasta la proxima arista que
 *  toque un borde incluyendo estas dos.
 */
double sumarBordes(AristaVoronoi **a, PuntoVoronoi *p1, PuntoVoronoi *p2, 
		   PuntoVoronoi *q1, PuntoVoronoi *q2, 
		   PoligonoVoronoi *pol) { 
  int tope1, tope2, i; 
  double suma = 0;
  double u1, v1, u2, v2;
  PuntoVoronoi *paux;

  if (DEBUG) {
    printf ("(q1x, q1y) = (%g, %g)\n", q1->x, q1->y);
    printf ("(q2x, q2y) = (%g, %g)\n", q2->x, q2->y);
  }
  
  if (p1->x > p2->x) { // bordeCuadrado los dio vuelta
    paux = q1;
    q1 = q2;
    q2 = paux;
  }
  u2 = q2->x;
  v2 = q2->y;
  
  tope1 = tope (q2);
  
  if (DEBUG) {
    printf("tope %d: (%g, %g)\n", tope1, u2, v2);
    printf("buscando otro tope\n");
  }
  do{	
    if ((*a)->poliDer == pol) {
      *a = (*a)->cwPred;
    }
    else if ((*a)->poliIzq == pol) {
      *a = (*a)->cwSucc;
    }
    if ((*a)->poliDer == pol) {
      p1 = (*a)->ptoFin;
      p2 = (*a)->ptoIni;
    }
    else if ((*a)->poliIzq == pol) {
      p1 = (*a)->ptoIni;
      p2 = (*a)->ptoFin;
    }
    if (DEBUG) {
      imprimirArista(*a);
    }
  } while(!bordeCuadrado (p1, p2, q1, q2));
  
  if (p1->x > p2->x) { // bordeCuadrado los dio vuelta
    paux = q1;
    q1 = q2;
    q2 = paux;
  }
  tope2 = tope (q1);
  
  if (DEBUG) {
    printf ("p1: (%g, %g)\n", p1->x, p1->y);
    printf ("p2: (%g, %g)\n", p2->x, p2->y);
    printf ("q1: (%g, %g)\n", q1->x, q1->y);
    printf ("q2: (%g, %g)\n", q2->x, q2->y);
    printf("encontrado otro tope %d: ", tope2);
    imprimirArista(*a);
  }
  for (i = tope1; ((i - 1) % 4 + 1) != tope2; i++) {
    u1 = u2;
    v1 = v2;
    switch ((i - 1) % 4 + 1) {
    case BORDE_INF:
      u2 = 1;
      v2 = 0;
      break;
    case BORDE_DER:
      u2 = 1;
      v2 = 1;
      break;
    case BORDE_SUP:
      u2 = 0;
      v2 = 1;
      break;
    case BORDE_IZQ:
      u2 = 0;
      v2 = 0;
      break;
    }
    suma += u1 * v2 - u2 * v1;
    if (DEBUG2) printf("sumando1 (%g, %g),(%g, %g)\n", u1, v1, u2, v2);
  }
  suma += u2 * q1->y - q1->x * v2;
  if (DEBUG2) printf("sumando2 (%g, %g),(%g, %g)\n", u2, v2, q1->x, q1->y);

  return suma;
}
