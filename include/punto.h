#ifndef VIS_PUNTO
#define VIS_PUNTO

#include "estructuras_malla.h"

#ifdef __cplusplus
extern "C" {
#endif  

struct puntoVoronoi {
  AristaVoronoi *a;
  double x, y;
  NodoLista *n;	
};

PuntoVoronoi *newPuntoVoronoi(double x, double y);

#ifdef __cplusplus
}
#endif  

#endif
