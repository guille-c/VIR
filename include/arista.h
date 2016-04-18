#ifndef VIS_ARISTA
#define VIS_ARISTA

#include "estructuras_malla.h"
#include "punto.h"
#include "poligono.h"

#ifdef __cplusplus
extern "C" {
#endif  

struct aristaVoronoi {
  PuntoVoronoi *ptoIni, *ptoFin;
  PoligonoVoronoi *poliDer, *poliIzq;
  AristaVoronoi *cwPred, *ccwPred, *cwSucc, *ccwSucc;
  NodoLista *n;
  double cosDer, cosIzq, sinDer, sinIzq; 
    // seno y coseno del angulo que forman los sitios derecho e
    // izquierdo con respecto al eje x.
};

AristaVoronoi *newAristaVoronoi(PuntoVoronoi *p1, 
				PuntoVoronoi *p2);
void eliminarAristaVoronoi(AristaVoronoi *a);
void imprimirArista(AristaVoronoi *a);
void calcularSinCos(AristaVoronoi *a);

#ifdef __cplusplus
}
#endif

#endif
