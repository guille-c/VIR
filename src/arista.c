#include "arista.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

AristaVoronoi *newAristaVoronoi(PuntoVoronoi *p1, 
				PuntoVoronoi *p2) {
  AristaVoronoi *a = 
    (AristaVoronoi * ) malloc(sizeof(AristaVoronoi));
  a->ptoIni = p1;
  a->ptoFin = p2;
  a->poliDer = NULL;
  a->poliIzq = NULL;
  a->cwPred = NULL;
  a->ccwPred = NULL;
  a->cwSucc = NULL;
  a->ccwSucc = NULL;
  a->cosDer = -2;
  a->cosIzq = -2;
  a->sinDer = -2;
  a->sinIzq = -2;
  return a;
}

void eliminarAristaVoronoi(AristaVoronoi *a) {
  free (a);
}

void imprimirArista(AristaVoronoi *a) {
  printf("{(%g, %g), (%g, %g)}\n",
	 a->ptoIni->x, a->ptoIni->y, a->ptoFin->x, a->ptoFin->y);
}

void calcularSinCos(AristaVoronoi *a) {
  PoligonoVoronoi *pol;
  double d;

  if (a == NULL || a->ptoIni == NULL || a->ptoFin == NULL) {
    fprintf (stderr, "ERROR en arista->calcularSinCos:");
    fprintf (stderr, " la arista o uno de sus ptos. no existe\n");
    exit (1);
  }
  
  d = sqrt (pow (a->ptoIni->x - a->ptoFin->x, 2) + 
	    pow(a->ptoIni->y - a->ptoFin->y, 2));

  a->sinDer = (a->ptoFin->x - a->ptoIni->x) / d;
  a->cosDer = (a->ptoIni->y - a->ptoFin->y) / d;

  a->sinIzq = (a->ptoIni->x - a->ptoFin->x) / d;
  a->cosIzq = (a->ptoFin->y - a->ptoIni->y) / d;
}

/*
void calcularSinCos(AristaVoronoi *a) {
  PoligonoVoronoi *pol;
  PuntoVoronoi *p1, *p2;
  double d;

  if (a == NULL || a->ptoIni == NULL || a->ptoFin == NULL) {
    fprintf (stderr, "ERROR en arista->calcularSinCos:");
    fprintf (stderr, " la arista o uno de sus ptos. no existe\n");
    exit (1);
  }
  
  d = sqrt (pow (a->ptoIni->x - a->ptoFin->x, 2) + 
	    pow(a->ptoIni->y - a->ptoFin->y, 2));

  if (a->poliDer == NULL) {
    a->cosDer = -2;
    a->sinDer = -2;
  }
  else {
    pol = a->poliDer;
    if ((a->ptoIni->x - pol->x) * (a->ptoFin->y - a->ptoIni->y) - 
	(a->ptoIni->y - pol->y) * (a->ptoFin->x - a->ptoIni->x) >= 0) {
      p1 = a->ptoIni;
      p2 = a->ptoFin;
    }
    else {
      p1 = a->ptoFin;
      p2 = a->ptoIni;
    }
    a->sinDer = (p1->x - p2->x) / d;
    a->cosDer = (p2->y - p1->y) / d;
  }

  if (a->poliIzq == NULL) {
    a->cosIzq = -2;
    a->sinIzq = -2;
  }
  else {
    pol = a->poliIzq;
    if ((a->ptoIni->x - pol->x) * (a->ptoFin->y - a->ptoIni->y) - 
	(a->ptoIni->y - pol->y) * (a->ptoFin->x - a->ptoIni->x) >= 0) {
      p1 = a->ptoIni;
      p2 = a->ptoFin;
    }
    else {
      p1 = a->ptoFin;
      p2 = a->ptoIni;
    }
    a->sinIzq = (p1->x - p2->x) / d;
    a->cosIzq = (p2->y - p1->y) / d;
  }
}
*/
