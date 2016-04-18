#include "punto.h"
#include <malloc.h>

PuntoVoronoi *newPuntoVoronoi(double x, double y) {
  PuntoVoronoi *p = (PuntoVoronoi *) malloc(sizeof(PuntoVoronoi));
  p->x = x;
  p->y = y;
  return p;
}
